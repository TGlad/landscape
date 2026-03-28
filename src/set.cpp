#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

static const double pi = std::acos(-1.0);

#include <iostream>
#include <vector>
#include <Eigen/Dense>

Landscape::Set::Set(const std::string &name, int num_balls) : name(name) 
{ 
  balls.resize(num_balls); 
  for (int i = 0; i<num_balls; i++)
    balls[i].parent_set = this;
  conn.resize(num_balls); // defaults to disconnected
} 
// multiple copy and asignment constructors that seem to be necessary
Landscape::Set::Set(const Set &o)
  : name(o.name), balls(o.balls), leaf_balls(o.leaf_balls),
    conn(o.conn), leaf_union(o.leaf_union)
{
  for (auto &b : balls)      b.parent_set = this;
  for (auto &b : leaf_balls) b.parent_set = this;
}
Landscape::Set::Set(Set &&o)
  : name(std::move(o.name)), balls(std::move(o.balls)),
    leaf_balls(std::move(o.leaf_balls)),
    conn(std::move(o.conn)), leaf_union(o.leaf_union)
{
  for (auto &b : balls)      b.parent_set = this;
  for (auto &b : leaf_balls) b.parent_set = this;
}
Landscape::Set &Landscape::Set::operator=(const Set &o)
{
  if (this == &o) return *this;
  name = o.name; balls = o.balls; leaf_balls = o.leaf_balls;
  conn = o.conn; leaf_union = o.leaf_union;
  for (auto &b : balls)      b.parent_set = this;
  for (auto &b : leaf_balls) b.parent_set = this;
  return *this;
}
Landscape::Set &Landscape::Set::operator=(Set &&o)
{
  if (this == &o) return *this;
  name = std::move(o.name); balls = std::move(o.balls);
  leaf_balls = std::move(o.leaf_balls);
  conn = std::move(o.conn); leaf_union = o.leaf_union;
  for (auto &b : balls)      b.parent_set = this;
  for (auto &b : leaf_balls) b.parent_set = this;
  return *this;
}

/* alternative method that seems to be simpler
void Landscape::Set::findOrthogonalSphere(int I, int J, int K, int L)
{
  Eigen::Matrix4d A;
  Eigen::Vector4d b;

  // We want a point P such that Power_0(P) = Power_1(P) = Power_2(P) = Power_3(P) = R^2
  // This forms a linear system in [Px, Py, Pz, -R^2]
  const std::array<int,4> idx = {I, J, K, L};
  for (int i = 0; i < 4; ++i) 
  {
    double k = balls[idx[i]].curvature;
    Eigen::Vector3d n = balls[idx[i]].dir;
    double d = balls[idx[i]].dist;

    double coeff_P = -2.0 * (k * d + 1.0);
    A(i, 0) = n.x() * coeff_P;
    A(i, 1) = n.y() * coeff_P;
    A(i, 2) = n.z() * coeff_P;
    A(i, 3) = k; // Coefficient for X

    b(i) = -(k * d * d + 2.0 * d);
  }

  // Solve for [Px, Py, Pz, X]
  Eigen::Vector4d sol = A.colPivHouseholderQr().solve(b);
  Eigen::Vector3d centre = sol.head<3>();
  double X = sol(3);
  
  // Since X = P^2 - R^2, then R^2 = P^2 - X
  double r_sq = centre.squaredNorm() - X;

  if (r_sq < 0.0)
    std::cerr << "Warning: orthogonal sphere is inside another sphere" << std::endl;

  float rad = std::sqrt(std::abs(r_sq));

  Ball leaf;
  leaf.parent_set = this;
  leaf.dir        = centre.normalized();
  leaf.dist       = centre.norm() - rad;
  leaf.curvature  = 1.0 / rad;
  leaf_balls.push_back(leaf);
}
*/

void Landscape::Set::addLeafBall(int i, int j, int k, int l)
{
  leaf_ball_ids.push_back(Eigen::Vector4i(i,j,k,l));
  leaf_ball_scales.push_back(1.0);
}
void Landscape::Set::addLeafBall(int i, int j, int k, double scale)
{
  leaf_ball_ids.push_back(Eigen::Vector4i(i,j,k,-1));
  leaf_ball_scales.push_back(scale);
}


void Landscape::Set::calculateLeafBall(int i, int j, int k, int l)
{
  // A sphere X orthogonal to sphere A satisfies |Cx-Ca|² = rx²+ra².
  // Expanding with w = |Cx|²-rx²:  2*Ca·Cx - w = |Ca|²-ra²
  // For a plane (curvature=0), 90° intersection means the center lies on the plane: n·Cx = dist
  // Four balls → 4×4 linear system in (Cx.x, Cx.y, Cx.z, w).

  const std::array<int,4> idx = {i, j, k, l};
  Eigen::Matrix4d M;
  Eigen::Vector4d rhs;

  for (int row = 0; row < 4; row++)
  {
    const Ball &b = balls[idx[row]];
    if (std::abs(b.curvature) > 1e-6)
    {
      double r = 1.0 / b.curvature;
      Eigen::Vector3d C = b.dir * (b.dist + r);
      M(row, 0) = 2.0 * C.x();
      M(row, 1) = 2.0 * C.y();
      M(row, 2) = 2.0 * C.z();
      M(row, 3) = -1.0;
      rhs(row) = C.squaredNorm() - r * r;
    }
    else
    {
      // Plane: center of orthogonal sphere lies on the plane → n·Cx = dist
      M(row, 0) = b.dir.x();
      M(row, 1) = b.dir.y();
      M(row, 2) = b.dir.z();
      M(row, 3) = 0.0;
      rhs(row) = b.dist;
    }
  }

  Eigen::Vector4d sol = M.fullPivLu().solve(rhs);
  Eigen::Vector3d Cx(sol(0), sol(1), sol(2));
  double w   = sol(3);           // w = |Cx|² - rx²
  double rx2 = Cx.squaredNorm() - w;

  if (rx2 <= 0.0) 
  {
    std::cerr << "[addLeafBall] degenerate: rx² = " << rx2 << " (no real orthogonal sphere)\n";
    return;
  }

  double rx = std::sqrt(rx2);
  Ball leaf;
  leaf.parent_set = this;
  leaf.dir        = Cx.normalized();
  if (leaf.dir == Eigen::Vector3d(0,0,0))
    leaf.dir[2] = 1.0;
  leaf.dist       = Cx.norm() - rx;
  leaf.curvature  = 1.0 / rx;
  leaf_balls.push_back(leaf);

  // ── Validation: check 90° intersection with each of the 4 balls ──────────
  bool ok = true;
  for (int row = 0; row < 4; row++)
  {
    const Ball &b = balls[idx[row]];
    double err;
    if (b.curvature != 0.0)
    {
      double rb = 1.0 / b.curvature;
      Eigen::Vector3d Cb = b.dir * (b.dist + rb);
      double d2 = (Cx - Cb).squaredNorm();
      // orthogonality: |Cx-Cb|² == rx²+rb²
      double target = rx2 + rb * rb;
      err = d2 - target;
    }
    else
    {
      // plane orthogonality: n·Cx == dist
      err = b.dir.dot(Cx) - b.dist;
    }
    bool ball_ok = std::abs(err) < 1e-6;
    if (!ball_ok)
      std::cout << "[addLeafBall] ball " << idx[row]
                << " ortho-err = " << err
                << (ball_ok ? "  OK" : "  FAIL") << "\n";
    if (!ball_ok) 
      ok = false;
  }
  if (!ok)
    std::cout << "[addLeafBall] rx = " << rx
              << "  |Cx| = " << Cx.norm()
              << "  " << (ok ? "ALL OK" : "FAIL") << "\n\n";
}

void Landscape::Set::calculateLeafBall(int i, int j, int k, double scale)
{
  // 4th row: Cx must lie in the plane of the 3 centres
  const std::array<int,3> idx = {i, j, k};
  std::vector<Eigen::Vector3d> centres(3);
  for (int row = 0; row < 3; row++)
  {
    const Ball &b = balls[idx[row]];
    if (b.curvature == 0.0)
      std::cerr << "Error: can't use this leafBall method on planes" << std::endl;
    centres[row] = b.dir * (b.dist + 1.0 / b.curvature);
  }

  Eigen::Vector3d v1 = centres[1] - centres[0];
  Eigen::Vector3d v2 = centres[2] - centres[0];
  Eigen::Vector3d n = v1.cross(v2);
  if (n.norm() < 1e-10)
  {
    std::cerr << "[addLeafBall3] degenerate: 3 centres are collinear\n";
    return;
  }
  double area = n.norm() / 2.0;
  double radius = std::sqrt(area / (2.0 * 3.14159));
  Eigen::Vector3d centre = centres[0] + (v2.squaredNorm()*n.cross(v1) + v1.squaredNorm()*v2.cross(n)) / (2.0 * n.squaredNorm());
  n.normalize();
  double dir = centre.dot(n) > 0.0 ? 1.0 : -1.0; // cheat which assumes set centred around 0,0,0
  centre -= n*radius*scale*dir;
  radius *= std::abs(scale);

  Ball leaf;
  leaf.parent_set = this;
  leaf.dir        = centre.normalized();
  leaf.dist       = centre.norm() - radius;
  leaf.curvature  = 1.0 / radius;
  leaf_balls.push_back(leaf);
}

bool Landscape::Set::verifyConnectivity(double tol) const
{
  bool all_pass = true;
  int n = (int)balls.size();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < i; j++)
    {
      int order = conn(i, j);
      const Ball &bi = balls[i];
      const Ball &bj = balls[j];

      double actual = 0, target = 0;
      const char *label = "angle";
      bool ok = false;

      if (order == 0)
      {
        // Separation: check gap >= k
        if (bi.curvature != 0.0 && bj.curvature != 0.0)
        {
          double ri = 1.0 / bi.curvature, rj = 1.0 / bj.curvature;
          Eigen::Vector3d Ci = bi.dir * (bi.dist + ri);
          Eigen::Vector3d Cj = bj.dir * (bj.dist + rj);
          actual = (Ci - Cj).norm();
          target = ri + rj; // k=0; gap = actual - target >= 0
          label = "gap";
          ok = (actual - target) >= -tol;
        }
        else if (bi.curvature != 0.0 || bj.curvature != 0.0)
        {
          const Ball &sphere = (bi.curvature != 0.0) ? bi : bj;
          const Ball &plane  = (bi.curvature != 0.0) ? bj : bi;
          double r = 1.0 / sphere.curvature;
          Eigen::Vector3d C = sphere.dir * (sphere.dist + r);
          actual = plane.dir.dot(C) - plane.dist;
          target = r; // sphere centre must be at least r beyond the plane
          label = "gap";
          ok = (actual - target) >= -tol;
        }
        else
        {
          continue; // plane-plane separation not meaningful
        }
        if (!ok)
        { 
          all_pass = false;
          std::cout << "  [" << name << "] balls (" << i << "," << j << ") order=0"
                    << " " << label << ": min=" << target << " actual=" << actual
                    << " gap=" << (actual - target) << (ok ? "  OK" : "  FAIL") << "\n";
        }
        continue;
      }

      if (bi.curvature != 0.0 && bj.curvature != 0.0)
      {
        double ri = 1.0 / bi.curvature, rj = 1.0 / bj.curvature;
        Eigen::Vector3d Ci = bi.dir * (bi.dist + ri);
        Eigen::Vector3d Cj = bj.dir * (bj.dist + rj);
        double d = (Ci - Cj).norm();
        if (order == -1)
        {
          label = "dist"; target = ri + rj; actual = d;
        }
        else
        {
          double cos_theta = (d*d - ri*ri - rj*rj) / (2.0 * ri * rj);
          target = pi / (double)order;
          actual = std::acos(std::clamp(cos_theta, -1.0, 1.0));
        }
      }
      else if (bi.curvature == 0.0 && bj.curvature == 0.0)
      {
        target = pi / (double)order;
        // Planes are treated as unoriented: opposing normals (n, -n) meet at angle 0,
        // parallel normals (n, n) meet at angle pi. Hence negate the dot product.
        actual = std::acos(std::clamp(-bi.dir.dot(bj.dir), -1.0, 1.0));
      }
      else
      {
        const Ball &sphere = (bi.curvature != 0.0) ? bi : bj;
        const Ball &plane  = (bi.curvature != 0.0) ? bj : bi;
        double r = 1.0 / sphere.curvature;
        Eigen::Vector3d C = sphere.dir * (sphere.dist + r);
        double cos_targ = (order == -1) ? 1.0 : std::cos(pi / (double)order);
        target = r * cos_targ;
        actual = plane.dir.dot(C) - plane.dist;
        label = "dist";
      }

      ok = std::abs(actual - target) <= tol;
      if (!ok) 
      {
        all_pass = false;
        std::cout << "  [" << name << "] balls (" << i << "," << j << ") order=" << order
                  << " " << label << ": target=" << target << " actual=" << actual
                  << " err=" << (actual - target) << (ok ? "  OK" : "  FAIL") << "\n";
      }
    }
  }
  return all_pass;
}