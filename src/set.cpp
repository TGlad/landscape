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
  // Expanding with w = |Cx|²-rx²:  2*Ca·Cx - w = |Ca|²-ra².
  // Four balls → 4×4 linear system in (Cx.x, Cx.y, Cx.z, w).

  const std::array<int,4> idx = {i, j, k, l};
  Eigen::Matrix4d M;
  Eigen::Vector4d rhs;

  for (int row = 0; row < 4; row++)
  {
    const Ball &b = balls[idx[row]];
    const double r = b.radius;
    const Eigen::Vector3d &C = b.centre;
    M(row, 0) = 2.0 * C.x();
    M(row, 1) = 2.0 * C.y();
    M(row, 2) = 2.0 * C.z();
    M(row, 3) = -1.0;
    rhs(row) = C.squaredNorm() - r * r;
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
  leaf.centre     = Cx;
  leaf.radius     = rx;
  leaf_balls.push_back(leaf);

  // ── Validation: check 90° intersection with each of the 4 balls ──────────
  bool ok = true;
  for (int row = 0; row < 4; row++)
  {
    const Ball &b = balls[idx[row]];
    double rb = b.radius;
    Eigen::Vector3d Cb = b.centre;
    double d2 = (Cx - Cb).squaredNorm();
    // orthogonality: |Cx-Cb|² == rx²+rb²
    double target = rx2 + rb * rb;
    double err = d2 - target;
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

bool triIntersectsSphere2(const std::vector<Eigen::Vector3d> &tri, const Eigen::Vector3d &centre, double rad, double eps)
{
  Eigen::Vector3d norm = (tri[2] - tri[0]).cross(tri[1]-tri[0]).normalized();
  if (tri[0].dot(norm) < 0.0)
    norm = -norm;
  double d_plane = tri[0].dot(norm);
  double d_centre = centre.dot(norm);
  return d_plane < d_centre + rad - eps;
}

bool triIntersectsSphere(const std::vector<Eigen::Vector3d> &tri, const Eigen::Vector3d &centre, double rad, double eps)
{
  if (tri.size() != 3 || rad <= 0.0)
    return false;

  const double margin = eps;
  const double effective_rad = rad - margin;
  if (effective_rad <= 0.0)
    return false;

  const Eigen::Vector3d &a = tri[0];
  const Eigen::Vector3d &b = tri[1];
  const Eigen::Vector3d &c = tri[2];

  // Closest point on triangle to sphere center (Ericson, Real-Time Collision Detection).
  Eigen::Vector3d ab = b - a;
  Eigen::Vector3d ac = c - a;
  Eigen::Vector3d ap = centre - a;

  double d1 = ab.dot(ap);
  double d2 = ac.dot(ap);
  if (d1 <= 0.0 && d2 <= 0.0)
    return (centre - a).squaredNorm() < effective_rad * effective_rad;

  Eigen::Vector3d bp = centre - b;
  double d3 = ab.dot(bp);
  double d4 = ac.dot(bp);
  if (d3 >= 0.0 && d4 <= d3)
    return (centre - b).squaredNorm() < effective_rad * effective_rad;

  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
    double v = d1 / (d1 - d3);
    Eigen::Vector3d closest = a + v * ab;
    return (centre - closest).squaredNorm() < effective_rad * effective_rad;
  }

  Eigen::Vector3d cp = centre - c;
  double d5 = ab.dot(cp);
  double d6 = ac.dot(cp);
  if (d6 >= 0.0 && d5 <= d6)
    return (centre - c).squaredNorm() < effective_rad * effective_rad;

  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
    double w = d2 / (d2 - d6);
    Eigen::Vector3d closest = a + w * ac;
    return (centre - closest).squaredNorm() < effective_rad * effective_rad;
  }

  double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
  {
    Eigen::Vector3d bc = c - b;
    double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    Eigen::Vector3d closest = b + w * bc;
    return (centre - closest).squaredNorm() < effective_rad * effective_rad;
  }

  // Inside face region.
  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom;
  double w = vc * denom;
  Eigen::Vector3d closest = a + ab * v + ac * w;
  return (centre - closest).squaredNorm() < effective_rad * effective_rad;
}

void Landscape::Set::calculateLeafBalls()
{
  // iterate leaf_ball_ids and run calculateLeafBall on each suitable quad
  // then filter out those outside the surface polyhedron

  for (int i = 0; i<(int)leaf_ball_set.size()-3; i++)
  {
    int I = leaf_ball_set[i];
    for (int j = i+1; j<(int)leaf_ball_set.size()-2; j++)
    {
      int J = leaf_ball_set[j];
      for (int k = j+1; k<(int)leaf_ball_set.size()-1; k++)
      {
        int K = leaf_ball_set[k];
        for (int l = k+1; l<(int)leaf_ball_set.size(); l++)
        {
          int L = leaf_ball_set[l];
          int count = conn(I,J)>0 ? 1 : 0;
          count += conn(I,K)>0 ? 1 : 0;
          count += conn(I,L)>0 ? 1 : 0;
          count += conn(J,K)>0 ? 1 : 0;
          count += conn(J,L)>0 ? 1 : 0;
          count += conn(K,L)>0 ? 1 : 0;
          if (count < 5)
            continue;
          calculateLeafBall(I,J,K,L);
        }
      }
    }
  }
  // OK now get the set of triangles:
  std::vector<Eigen::Vector3i> tris;
  double max_dist = 0.0;
  for (int i = 0; i<(int)leaf_ball_set.size()-2; i++)
  {
    int I = leaf_ball_set[i];
    Eigen::Vector3i tri;
    tri[0] = I;
    for (int j = i+1; j<(int)leaf_ball_set.size()-1; j++)
    {
      int J = leaf_ball_set[j];
      tri[1] = J;
      for (int k = j+1; k<(int)leaf_ball_set.size(); k++)
      {
        int K = leaf_ball_set[k];
        tri[2] = K;
        if (conn(I,J)>0 && conn(I,K)>0 && conn(J,K)>0)
        {
          for (int c = 0; c<3; c++)
          {
            const Ball &b = balls[tri[c]];
            max_dist = std::max(max_dist, b.centre.norm() + b.radius);  
          }
          tris.push_back(tri);
        }
      }
    }
  }
  // now if any triangle intersects any leaf ball then destroy the leaf ball and move on
  for (auto &tri: tris)
  {
    std::vector<Eigen::Vector3d> vs(3);
    for (int i = 0; i<3; i++)
    {
      const Ball &b = balls[tri[i]];
      vs[i] = b.centre;  
    }
    for (int i = leaf_balls.size()-1; i>=0; i--)
    {
      const Ball &b = leaf_balls[i];
      double r = b.radius;
      Eigen::Vector3d c = b.centre;
      float eps = 0.0001;
      if (triIntersectsSphere2(vs, c, r, eps) || c.norm() > max_dist) // max dist is a hack, need better option
      {
        leaf_balls[i] = leaf_balls.back();
        leaf_balls.pop_back();
      }
    }
  }
  for (auto &tri: tris)
  {
    if (conn(tri[0], tri[1]) != 3 || conn(tri[1],tri[2]) != 3 || conn(tri[0],tri[2]) != 3)
      continue;
    std::vector<Eigen::Vector3d> vs(3);
    for (int i = 0; i<3; i++)
    {
      const Ball &b = balls[tri[i]];
      vs[i] = b.centre;  
    }
    bool any_intersect_plane = false;
    for (int i = leaf_balls.size()-1; i>=0; i--)
    {
      const Ball &b = leaf_balls[i];
      double r = b.radius;
      Eigen::Vector3d c = b.centre;
      double eps = 0.0001;
      if (triIntersectsSphere2(vs, c, r, -eps))
        any_intersect_plane = true;
    }
    if (!any_intersect_plane)
    {
      calculateLeafBall(tri[0],tri[1],tri[2],0.5);      
    }
  }
}

Eigen::Vector3d findMeetingPoint(const std::vector<Eigen::Vector3d>& cs, const std::vector<double>& rs) 
{
// 1. Shift everything so cs[0] is the origin to improve numerical stability
  Eigen::Vector3d origin = cs[0];
  Eigen::Matrix<double, 2, 3> A;
  Eigen::Vector2d b;

  for (int i = 1; i < 3; ++i) 
  {
    Eigen::Vector3d relative_c = cs[i] - origin;
    // Linear equation: 2 * P . relative_c = r0^2 - ri^2 + ||relative_c||^2
    A.row(i - 1) = 2.0 * relative_c.transpose();  
    double r0_sq = rs[0] * rs[0];
    double ri_sq = rs[i] * rs[i];
    double dist_sq = relative_c.squaredNorm();
    b(i - 1) = r0_sq - ri_sq + dist_sq;
  }

  // 2. Solve the underdetermined system A * P = b
  // CompleteOrthogonalDecomposition finds the point P with the minimum norm.
  // Since our system is centered at cs[0] and the rows of A span the plane 
  // of the triangle, this solution is guaranteed to lie in that plane.
  Eigen::Vector3d p_relative = A.completeOrthogonalDecomposition().solve(b);

  // 3. Shift back to world coordinates
  return p_relative + origin;
}

void Landscape::Set::calculateLeafBall(int i, int j, int k, double scale)
{
  // 4th row: Cx must lie in the plane of the 3 centres
  const std::array<int,3> idx = {i, j, k};
  std::vector<Eigen::Vector3d> cs(3);
  std::vector<double> rs(3);
  for (int row = 0; row < 3; row++)
  {
    const Ball &b = balls[idx[row]];
    cs[row] = b.centre;
    rs[row] = b.radius;
  }

  Eigen::Vector3d v1 = cs[1] - cs[0];
  Eigen::Vector3d v2 = cs[2] - cs[0];
  Eigen::Vector3d n = v1.cross(v2);
  double area = n.norm() / 2.0;
  n.normalize();
  double radius = std::sqrt(area / (2.0 * 3.14159));

  Eigen::Vector3d centre = findMeetingPoint(cs, rs);
  double dir = centre.dot(n) > 0.0 ? 1.0 : -1.0; // cheat which assumes set centred around 0,0,0
  centre -= n*radius*scale*dir;
  radius *= std::abs(scale);

  Ball leaf;
  leaf.parent_set = this;
  leaf.centre     = centre;
  leaf.radius     = radius;
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
        double ri = bi.radius, rj = bj.radius;
        actual = (bi.centre - bj.centre).norm();
        target = ri + rj; // k=0; gap = actual - target >= 0
        label = "gap";
        ok = (actual - target) >= -tol;
        if (!ok)
        { 
          all_pass = false;
          std::cout << "  [" << name << "] balls (" << i << "," << j << ") order=0"
                    << " " << label << ": min=" << target << " actual=" << actual
                    << " gap=" << (actual - target) << (ok ? "  OK" : "  FAIL") << "\n";
        }
        continue;
      }

      double ri = bi.radius, rj = bj.radius;
      double d = (bi.centre - bj.centre).norm();
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