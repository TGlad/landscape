#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

static const double pi = std::acos(-1.0);

void Landscape::Set::addLeafBall(int i, int j, int k, int l)
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
    if (b.curvature != 0.0)
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
    std::cout << "[addLeafBall] ball " << idx[row]
              << " ortho-err = " << err
              << (ball_ok ? "  OK" : "  FAIL") << "\n";
    if (!ball_ok) ok = false;
  }
  std::cout << "[addLeafBall] rx = " << rx
            << "  |Cx| = " << Cx.norm()
            << "  " << (ok ? "ALL OK" : "FAIL") << "\n\n";
}

void Landscape::Set::addLeafBall(int i, int j, int k)
{
  // Orthogonal sphere whose centre lies in the plane of the 3 ball centres.
  // 3 orthogonality rows (same form as the 4-ball version) + 1 plane-constraint row.

  const std::array<int,3> idx = {i, j, k};
  Eigen::Matrix4d M;
  Eigen::Vector4d rhs;

  std::array<Eigen::Vector3d, 3> centers;
  for (int row = 0; row < 3; row++)
  {
    const Ball &b = balls[idx[row]];
    if (b.curvature != 0.0)
    {
      double r = 1.0 / b.curvature;
      Eigen::Vector3d C = b.dir * (b.dist + r);
      centers[row] = C;
      M(row, 0) = 2.0 * C.x();
      M(row, 1) = 2.0 * C.y();
      M(row, 2) = 2.0 * C.z();
      M(row, 3) = -1.0;
      rhs(row) = C.squaredNorm() - r * r;
    }
    else
    {
      // Plane: orthogonal sphere centre lies on the plane; use closest point to origin for centre
      centers[row] = b.dir * b.dist;
      M(row, 0) = b.dir.x();
      M(row, 1) = b.dir.y();
      M(row, 2) = b.dir.z();
      M(row, 3) = 0.0;
      rhs(row) = b.dist;
    }
  }

  // 4th row: Cx must lie in the plane of the 3 centres
  Eigen::Vector3d n = (centers[1] - centers[0]).cross(centers[2] - centers[0]);
  if (n.norm() < 1e-10)
  {
    std::cerr << "[addLeafBall3] degenerate: 3 centres are collinear\n";
    return;
  }
  n.normalize();
  double plane_d = n.dot(centers[0]);
  M(3, 0) = n.x();
  M(3, 1) = n.y();
  M(3, 2) = n.z();
  M(3, 3) = 0.0;
  rhs(3) = plane_d;

  Eigen::Vector4d sol = M.fullPivLu().solve(rhs);
  Eigen::Vector3d Cx(sol(0), sol(1), sol(2));
  double w   = sol(3);
  double rx2 = Cx.squaredNorm() - w;

  if (rx2 <= 0.0)
  {
    std::cerr << "[addLeafBall3] degenerate: rx\u00b2 = " << rx2 << " (no real orthogonal sphere)\n";
    return;
  }

  double rx = std::sqrt(rx2);
  Ball leaf;
  leaf.parent_set = this;
  leaf.dir        = Cx.normalized();
  leaf.dist       = Cx.norm() - rx;
  leaf.curvature  = 1.0 / rx;
  leaf_balls.push_back(leaf);
}

void Landscape::Set::applyConnectivity()
{
  // Gauss-Seidel iterated least squares: for each pair constraint compute the gradient of the
  // error w.r.t. the state (dir, dist, curvature) and apply the minimum-norm correction.
  // Kissing (order=-1) uses a distance constraint to avoid 1/sin(0) instability.
  // Sphere-plane uses a linear signed-distance constraint (always well-posed).
  const double damping = 1e-10;
  const double k = 0.0; // minimum required gap beyond ri+rj for disjoint pairs

  // Build a list of all active constraint pairs and shuffle each iteration
  // to avoid Gauss-Seidel ordering bias (earlier-indexed balls dominating).
  int n = (int)balls.size();
  std::vector<std::pair<int,int>> pairs;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      pairs.push_back({i, j});
  std::mt19937 rng(42);

  for (int it = 0; it < 2000; it++)
  {
    std::shuffle(pairs.begin(), pairs.end(), rng);
    for (auto [i, j] : pairs)
    {
        int order = conn(i, j);

        Ball &bi = balls[i];
        Ball &bj = balls[j];

        // Gradients of error w.r.t. state: de/d(dir) [tangent], de/d(dist), de/d(curvature)
        Eigen::Vector3d g_dir_i = Eigen::Vector3d::Zero(), g_dir_j = Eigen::Vector3d::Zero();
        double g_dist_i = 0, g_curv_i = 0, g_dist_j = 0, g_curv_j = 0;
        double error = 0;

        if (bi.curvature != 0.0 && bj.curvature != 0.0)
        {
          // ── Sphere–sphere ──────────────────────────────────────────────────
          double ri = 1.0 / bi.curvature, rj = 1.0 / bj.curvature;
          Eigen::Vector3d Ci = bi.dir * (bi.dist + ri);
          Eigen::Vector3d Cj = bj.dir * (bj.dist + rj);
          Eigen::Vector3d Delta = Ci - Cj;
          double d = Delta.norm();
          if (d < 1e-12) continue;

          if (order <= 0)
          {
            // Distance constraint: covers both kissing (order=-1, targ=ri+rj)
            // and separation (order=0, targ=ri+rj+k, unilateral).
            double targ_d = ri + rj + (order == 0 ? k : 0.0);
            error = d - targ_d;
            if (order == 0 && error >= 0.0) continue; // unilateral: skip if already separated
            Eigen::Vector3d dddCi =  Delta / d;
            Eigen::Vector3d dddCj = -Delta / d;
            g_dist_i = dddCi.dot(bi.dir);
            g_dist_j = dddCj.dot(bj.dir);
            // de/dkappa: chain through both C (via r=1/kappa) and the target
            g_curv_i = (dddCi.dot(bi.dir) - 1.0) * (-1.0 / (bi.curvature * bi.curvature));
            g_curv_j = (dddCj.dot(bj.dir) - 1.0) * (-1.0 / (bj.curvature * bj.curvature));
            g_dir_i = (bi.dist + ri) * (dddCi - dddCi.dot(bi.dir) * bi.dir);
            g_dir_j = (bj.dist + rj) * (dddCj - dddCj.dot(bj.dir) * bj.dir);
          }
          else
          {
            // Angle constraint: θ = π/order
            // cos θ = f = (d² - ri² - rj²) / (2 ri rj)
            double d2 = d * d;
            double cos_theta = (d2 - ri*ri - rj*rj) / (2.0 * ri * rj);
            if (cos_theta >= 1.0)
            {
              // Spheres don't intersect yet (too far apart): pull them to kissing distance
              // as a fallback so subsequent iterations can apply the angle constraint.
              error = d - (ri + rj);
              Eigen::Vector3d dddCi =  Delta / d;
              Eigen::Vector3d dddCj = -Delta / d;
              g_dist_i = dddCi.dot(bi.dir);
              g_dist_j = dddCj.dot(bj.dir);
              g_curv_i = (dddCi.dot(bi.dir) - 1.0) * (-1.0 / (bi.curvature * bi.curvature));
              g_curv_j = (dddCj.dot(bj.dir) - 1.0) * (-1.0 / (bj.curvature * bj.curvature));
              g_dir_i = (bi.dist + ri) * (dddCi - dddCi.dot(bi.dir) * bi.dir);
              g_dir_j = (bj.dist + rj) * (dddCj - dddCj.dot(bj.dir) * bj.dir);
            }
            else if (cos_theta <= -1.0)
            {
              // One sphere inside the other: push them apart to kissing distance
              error = d - (ri + rj);
              Eigen::Vector3d dddCi =  Delta / d;
              Eigen::Vector3d dddCj = -Delta / d;
              g_dist_i = dddCi.dot(bi.dir);
              g_dist_j = dddCj.dot(bj.dir);
              g_curv_i = (dddCi.dot(bi.dir) - 1.0) * (-1.0 / (bi.curvature * bi.curvature));
              g_curv_j = (dddCj.dot(bj.dir) - 1.0) * (-1.0 / (bj.curvature * bj.curvature));
              g_dir_i = (bi.dist + ri) * (dddCi - dddCi.dot(bi.dir) * bi.dir);
              g_dir_j = (bj.dist + rj) * (dddCj - dddCj.dot(bj.dir) * bj.dir);
            }
            else
            {
              double theta = std::acos(cos_theta);
              double sin_theta = std::sin(theta);
              if (std::abs(sin_theta) < 1e-10) continue;
              error = pi / (double)order - theta;
              double inv_sin = 1.0 / sin_theta;
              Eigen::Vector3d dfdCi =  Delta / (ri * rj);  // df/dCi = Δ/(ri rj)
              Eigen::Vector3d dfdCj = -dfdCi;
              double dfdri = -(ri*ri + d2 - rj*rj) / (2.0 * ri*ri * rj);
              double dfdrj = -(rj*rj + d2 - ri*ri) / (2.0 * rj*rj * ri);
              g_dist_i = inv_sin * dfdCi.dot(bi.dir);
              g_dist_j = inv_sin * dfdCj.dot(bj.dir);
              g_curv_i = inv_sin * (dfdri + dfdCi.dot(bi.dir)) * (-1.0 / (bi.curvature * bi.curvature));
              g_curv_j = inv_sin * (dfdrj + dfdCj.dot(bj.dir)) * (-1.0 / (bj.curvature * bj.curvature));
              g_dir_i = inv_sin * (bi.dist + ri) * (dfdCi - dfdCi.dot(bi.dir) * bi.dir);
              g_dir_j = inv_sin * (bj.dist + rj) * (dfdCj - dfdCj.dot(bj.dir) * bj.dir);
            } // end intersecting else
          }
        }
        else if (bi.curvature == 0.0 && bj.curvature == 0.0)
        {
          // ── Plane–plane ────────────────────────────────────────────────────
          // Kissing parallel planes and separation are degenerate in 3D; only angle makes sense.
          if (order <= 0) continue;
          double dot = bi.dir.dot(bj.dir);
          double theta = std::acos(std::clamp(dot, -1.0, 1.0));
          double sin_theta = std::sin(theta);
          if (std::abs(sin_theta) < 1e-10) continue;
          error = pi / (double)order - theta;
          double inv_sin = 1.0 / sin_theta;
          g_dir_i = inv_sin * (bj.dir - dot * bi.dir);
          g_dir_j = inv_sin * (bi.dir - dot * bj.dir);
          // dist and curvature don't enter the plane-plane angle
        }
        else
        {
          // ── Sphere–plane ───────────────────────────────────────────────────
          // All orders use the same signed-distance constraint:
          //   n·C - d_plane = r·cos_targ
          // cos_targ encodes the target angle (1 for kissing/separation, cos(π/order) otherwise).
          // For separation (order=0) an extra gap k is added and the constraint is unilateral.
          const bool i_is_sphere = (bi.curvature != 0.0);
          Ball &sphere = i_is_sphere ? bi : bj;
          Ball &plane  = i_is_sphere ? bj : bi;
          double r = 1.0 / sphere.curvature;
          Eigen::Vector3d C = sphere.dir * (sphere.dist + r);
          double signed_dist = plane.dir.dot(C) - plane.dist;
          double cos_targ = (order <= 0) ? 1.0 : std::cos(pi / (double)order);
          double targ_dist = r * cos_targ + (order == 0 ? k : 0.0);
          error = signed_dist - targ_dist;
          if (order == 0 && error >= 0.0) continue; // unilateral separation

          double g_dist_s = plane.dir.dot(sphere.dir);
          // d(targ_dist)/dr = cos_targ (the +k term is independent of r)
          double g_curv_s = (plane.dir.dot(sphere.dir) - cos_targ) * (-1.0 / (sphere.curvature * sphere.curvature));
          Eigen::Vector3d g_dir_s = (sphere.dist + r) * (plane.dir - plane.dir.dot(sphere.dir) * sphere.dir);
          double g_dist_p = -1.0;
          Eigen::Vector3d g_dir_p = C - plane.dir.dot(C) * plane.dir;

          if (i_is_sphere) 
          {
            g_dist_i = g_dist_s; g_curv_i = g_curv_s; g_dir_i = g_dir_s;
            g_dist_j = g_dist_p; g_curv_j = 0.0;      g_dir_j = g_dir_p;
          } 
          else 
          {
            g_dist_j = g_dist_s; g_curv_j = g_curv_s; g_dir_j = g_dir_s;
            g_dist_i = g_dist_p; g_curv_i = 0.0;      g_dir_i = g_dir_p;
          }
        }

        // Minimum-norm correction: δstate = -(error / ||g||²) * g
        double g2 = g_dir_i.squaredNorm() + g_dist_i*g_dist_i + g_curv_i*g_curv_i
                  + g_dir_j.squaredNorm() + g_dist_j*g_dist_j + g_curv_j*g_curv_j;
        double step = -error / (g2 + damping);

        bi.dir       += step * g_dir_i;  bi.dir.normalize();
        bi.dist      += step * g_dist_i;
        bi.curvature += step * g_curv_i;
        bj.dir       += step * g_dir_j;  bj.dir.normalize();
        bj.dist      += step * g_dist_j;
        bj.curvature += step * g_curv_j;
    } // end pairs loop
  } // end iterations
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
        if (!ok) all_pass = false;
        std::cout << "  [" << name << "] balls (" << i << "," << j << ") order=0"
                  << " " << label << ": min=" << target << " actual=" << actual
                  << " gap=" << (actual - target) << (ok ? "  OK" : "  FAIL") << "\n";
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
        actual = std::acos(std::clamp(bi.dir.dot(bj.dir), -1.0, 1.0));
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
      if (!ok) all_pass = false;
      std::cout << "  [" << name << "] balls (" << i << "," << j << ") order=" << order
                << " " << label << ": target=" << target << " actual=" << actual
                << " err=" << (actual - target) << (ok ? "  OK" : "  FAIL") << "\n";
    }
  }
  return all_pass;
}

void Landscape::outputCode(const std::string &filename) const
{
  std::ofstream out(filename);
  if (!out) { std::cerr << "Failed to open " << filename << "\n"; return; }

  // Sanitize a set name into a valid GLSL identifier
  auto ident = [](const std::string &name) {
    std::string s = name;
    for (char &c : s) if (!std::isalnum(c)) c = '_';
    return s;
  };

  // Build per-set offsets into the flat ball array
  int num_sets = (int)sets.size();
  std::vector<int> offsets(num_sets + 1, 0);
  for (int si = 0; si < num_sets; si++)
    offsets[si + 1] = offsets[si] + (int)sets[si].balls.size();
  int total_balls = offsets[num_sets];

  // Build per-set offsets into the flat leaf_ball array
  std::vector<int> leaf_offsets(num_sets + 1, 0);
  for (int si = 0; si < num_sets; si++)
    leaf_offsets[si + 1] = leaf_offsets[si] + (int)sets[si].leaf_balls.size();
  int total_leaf_balls = leaf_offsets[num_sets];

  // Find the flat index of a Ball pointer
  auto flatIndex = [&](const Set::Ball *ball) -> int {
    for (int si = 0; si < num_sets; si++)
      for (int bi = 0; bi < (int)sets[si].balls.size(); bi++)
        if (&sets[si].balls[bi] == ball) return offsets[si] + bi;
    return -1;
  };

  out << std::fixed << std::setprecision(7);

  // ── Ball struct ──────────────────────────────────────────────────────────
  out << "struct Ball {\n"
      << "    vec3  dir;\n"
      << "    float dist;\n"
      << "    float curvature;\n"
      << "    int   dest_ball;  // flat index into BALLS[]; -1 = none\n"
      << "};\n\n";

  // ── Set offset table ─────────────────────────────────────────────────────
  // SET_OFFSET[i] = first index in BALLS[] for set i
  // Ball b = BALLS[SET_OFFSET[dest_set] + local_ball_index]
  out << "const int NUM_SETS = " << num_sets << ";\n";
  out << "const int SET_OFFSET[" << num_sets << "] = int[" << num_sets << "](";
  for (int si = 0; si < num_sets; si++)
    out << offsets[si] << (si < num_sets - 1 ? ", " : "");
  out << ");\n";
  out << "const int SET_SIZE[" << num_sets << "] = int[" << num_sets << "](";
  for (int si = 0; si < num_sets; si++)
    out << (int)sets[si].balls.size() << (si < num_sets - 1 ? ", " : "");
  out << ");\n";
  out << "const int LEAF_OFFSET[" << num_sets << "] = int[" << num_sets << "](";
  for (int si = 0; si < num_sets; si++)
    out << leaf_offsets[si] << (si < num_sets - 1 ? ", " : "");
  out << ");\n";
  out << "const int LEAF_SIZE[" << num_sets << "] = int[" << num_sets << "](";
  for (int si = 0; si < num_sets; si++)
    out << (int)sets[si].leaf_balls.size() << (si < num_sets - 1 ? ", " : "");
  out << ");\n";
  out << "const bool LEAF_UNION[" << num_sets << "] = bool[" << num_sets << "](";
  for (int si = 0; si < num_sets; si++)
    out << (sets[si].leaf_union ? "true" : "false") << (si < num_sets - 1 ? ", " : "");
  out << ");\n\n";

  // ── Flat ball array ───────────────────────────────────────────────────────
  out << "const int NUM_BALLS = " << total_balls << ";\n";
  out << "const Ball BALLS[" << total_balls << "] = Ball[" << total_balls << "](\n";
  for (int si = 0; si < num_sets; si++)
  {
    const Set &s = sets[si];
    int n = (int)s.balls.size();
    out << "    // set " << si << ": " << s.name << "\n";
    for (int bi = 0; bi < n; bi++)
    {
      const Set::Ball &b = s.balls[bi];
      int dest = flatIndex(b.dest_ball);
      bool last = (si == num_sets - 1 && bi == n - 1);
      out << "    Ball(vec3("
          << b.dir.x() << ", " << b.dir.y() << ", " << b.dir.z() << "), "
          << b.dist << ", "
          << b.curvature << ", "
          << dest << ")"
          << (last ? "" : ",") << "\n";
    }
  }
  out << ");\n\n";

  // ── Flat leaf_ball array ──────────────────────────────────────────────────
  if (total_leaf_balls > 0)
  {
    out << "const int NUM_LEAF_BALLS = " << total_leaf_balls << ";\n";
    out << "const Ball LEAF_BALLS[" << total_leaf_balls << "] = Ball[" << total_leaf_balls << "](\n";
    int flat = 0;
    for (int si = 0; si < num_sets; si++)
    {
      const Set &s = sets[si];
      int n = (int)s.leaf_balls.size();
      if (n == 0) continue;
      out << "    // set " << si << ": " << s.name << "\n";
      for (int bi = 0; bi < n; bi++, flat++)
      {
        const Set::Ball &b = s.leaf_balls[bi];
        out << "    Ball(vec3("
            << b.dir.x() << ", " << b.dir.y() << ", " << b.dir.z() << "), "
            << b.dist << ", "
            << b.curvature << ", "
            << "-1)"  // leaf balls don't recurse
            << (flat < total_leaf_balls - 1 ? "," : "") << "\n";
      }
    }
    out << ");\n\n";
  }

  // ── Per-set connectivity ──────────────────────────────────────────────────
  for (int si = 0; si < num_sets; si++)
  {
    const Set &s = sets[si];
    int n = (int)s.balls.size();
    int nadj = n * (n + 1) / 2;
    std::string id = ident(s.name);

    out << "// ── Set " << si << ": " << s.name << " ────────────────────────────────────\n";
    out << "const int CONN_" << id << "[" << nadj << "] = int[" << nadj << "](";
    for (int k = 0; k < nadj; k++)
      out << s.conn.data[k] << (k < nadj - 1 ? ", " : "");
    out << ");\n";

    out << "int conn_" << id << "(int i, int j) {\n"
        << "    return i >= j ? CONN_" << id << "[i*(i+1)/2 + j]\n"
        << "                  : CONN_" << id << "[j*(j+1)/2 + i];\n"
        << "}\n\n";
  }

  // ── Möbius transforms ─────────────────────────────────────────────────────
  // Only balls with a dest_set carry a non-identity transform.  We build a
  // compact array of those transforms and an indirection table BALL_MOBIUS[]
  // that maps a flat ball index to its Möbius index, or -1 if it has none.
  //
  // GLSL usage:
  //   int mi = BALL_MOBIUS[ball_idx];
  //   vec3 p2 = (mi >= 0) ? applyMobius(mi, p) : p;

  // Collect flat indices of balls that need a transform.
  std::vector<int> mobius_flat; // flat ball indices that have dest_set
  std::vector<int> ball_mobius(total_balls, -1); // flat index → mobius index
  for (int si = 0; si < num_sets; si++) {
    const Set &s = sets[si];
    for (int bi = 0; bi < (int)s.balls.size(); bi++) {
      if (!s.balls[bi].dest_set.empty()) {
        ball_mobius[offsets[si] + bi] = (int)mobius_flat.size();
        mobius_flat.push_back(offsets[si] + bi);
      }
    }
  }
  int num_mobius = (int)mobius_flat.size();

  // BALL_MOBIUS indirection table
  out << "const int BALL_MOBIUS[" << total_balls << "] = int[" << total_balls << "](";
  for (int i = 0; i < total_balls; i++)
    out << ball_mobius[i] << (i < total_balls-1 ? "," : "");
  out << ");\n\n";

  if (num_mobius > 0)
  {
    // Helper to look up (set index, ball index) from a flat index
    auto ballAt = [&](int flat) -> const Set::Ball & {
      for (int si = 0; si < num_sets; si++)
        if (flat < offsets[si+1])
          return sets[si].balls[flat - offsets[si]];
      return sets[0].balls[0]; // unreachable
    };

    // MOBIUS_C — inversion centres (unused when MOBIUS_SIM is true)
    out << "const vec3 MOBIUS_C[" << num_mobius << "] = vec3[" << num_mobius << "](\n";
    for (int mi = 0; mi < num_mobius; mi++) {
      const auto &C = ballAt(mobius_flat[mi]).mobius.C;
      out << "    vec3(" << C.x() << "," << C.y() << "," << C.z() << ")"
          << (mi < num_mobius-1 ? "," : "") << "\n";
    }
    out << ");\n\n";

    // MOBIUS_T — translations / images of ∞
    out << "const vec3 MOBIUS_T[" << num_mobius << "] = vec3[" << num_mobius << "](\n";
    for (int mi = 0; mi < num_mobius; mi++) {
      const auto &T = ballAt(mobius_flat[mi]).mobius.T;
      out << "    vec3(" << T.x() << "," << T.y() << "," << T.z() << ")"
          << (mi < num_mobius-1 ? "," : "") << "\n";
    }
    out << ");\n\n";

    // MOBIUS_S — scale factors
    out << "const float MOBIUS_S[" << num_mobius << "] = float[" << num_mobius << "](";
    for (int mi = 0; mi < num_mobius; mi++)
      out << ballAt(mobius_flat[mi]).mobius.s << (mi < num_mobius-1 ? "," : "");
    out << ");\n\n";

    // MOBIUS_R — rotation/reflection matrices (GLSL mat3 is column-major)
    out << "const mat3 MOBIUS_R[" << num_mobius << "] = mat3[" << num_mobius << "](\n";
    for (int mi = 0; mi < num_mobius; mi++) {
      const Eigen::Matrix3d &Rm = ballAt(mobius_flat[mi]).mobius.R;
      out << "    mat3("
          << Rm(0,0) << "," << Rm(1,0) << "," << Rm(2,0) << ","
          << Rm(0,1) << "," << Rm(1,1) << "," << Rm(2,1) << ","
          << Rm(0,2) << "," << Rm(1,2) << "," << Rm(2,2) << ")"
          << (mi < num_mobius-1 ? "," : "") << "\n";
    }
    out << ");\n\n";

    // MOBIUS_SIM — true iff the transform is a pure similarity (no inversion)
    out << "const bool MOBIUS_SIM[" << num_mobius << "] = bool[" << num_mobius << "](";
    for (int mi = 0; mi < num_mobius; mi++)
      out << (ballAt(mobius_flat[mi]).mobius.is_similarity ? "true" : "false")
          << (mi < num_mobius-1 ? "," : "");
    out << ");\n\n";

    // applyMobius() GLSL function
    out << "// Apply the Möbius transform for Möbius index mi to point v.\n"
        << "//   Similarity (MOBIUS_SIM=true):  v' = T + s * R * v\n"
        << "//   Inversion  (MOBIUS_SIM=false): v' = T + (s/dot(w,w)) * R * w,  w = v - C\n"
        << "vec3 applyMobius(int mi, vec3 v) {\n"
        << "    if (MOBIUS_SIM[mi])\n"
        << "        return MOBIUS_T[mi] + MOBIUS_S[mi] * (MOBIUS_R[mi] * v);\n"
        << "    vec3 w = v - MOBIUS_C[mi];\n"
        << "    return MOBIUS_T[mi] + (MOBIUS_S[mi] / dot(w, w)) * (MOBIUS_R[mi] * w);\n"
        << "}\n\n";
  }

  std::cout << "Wrote " << filename << "\n";
}

void Landscape::addSetToTypes(Set &set)
{
  int set_size = (int)set.balls.size();
  for (int i = 0; i < set_size; i++)
  {
    std::vector<int> out_is;
    Type new_type(set.conn, i, out_is);

    // Populate per-ball permutation maps.
    Set::Ball &ball = set.balls[i];
    ball.type_to_set = out_is;
    ball.set_to_type.assign(set_size, -1);
    for (int ti = 0; ti < (int)out_is.size(); ti++)
      ball.set_to_type[out_is[ti]] = ti;

    bool found = false;
    for (int j = 0; j < (int)types.size(); j++)
    {
      if (types[j].conn.data == new_type.conn.data)
      {
        // new_type already exists, so add it in.
        types[j].balls.push_back(&ball);
        ball.type = &types[j];
        found = true;
        break;
      }
    }
    if (!found)
    {
      types.push_back(new_type);
      types.back().balls.push_back(&ball);
      ball.type = &types.back();
    }
  }
}

// ── Conformal (O(4,1)) Möbius transform ──────────────────────────────────────
//
// Metric η = diag(1,1,1,1,-1) on ℝ^{4,1}.
// Sphere (c,r) → σ̂ = (c, (1-|c|²+r²)/2, (1+|c|²-r²)/2) / r   (σ̂·σ̂ = 1).
// Point  p    → P  = (p, (1-|p|²)/2, (1+|p|²)/2)              (P·P  = 0).
// M ∈ O(4,1) satisfies M^T η M = η and maps source σ̂s to dest σ̂s.
// With exactly 5 independent sphere pairs M = C5 * A5^{-1} is exact.

using Vec5 = Eigen::Matrix<double,5,1>;
using Mat5 = Eigen::Matrix<double,5,5>;

static double mink_dot(const Vec5 &a, const Vec5 &b)
{
  return a.head<4>().dot(b.head<4>()) - a(4)*b(4);
}

static Vec5 conformal_sphere(const Eigen::Vector3d &c, double r)
{
  double c2 = c.squaredNorm();
  Vec5 v;
  v << c.x(), c.y(), c.z(), (1.0-c2+r*r)/2.0, (1.0+c2-r*r)/2.0;
  return v / r; // unit spacelike: v·η·v = 1
}

// Decompose an O(4,1) matrix M into the GLSL-friendly T,C,s,R form stored in
// Landscape::Set::Ball::Mobius.
//
// Every element of O(4,1) acts on ℝ³∪{∞} as either:
//   Similarity  (∞ → ∞):   f(v) = T + s·R·v                    [is_similarity=true]
//   Inversion   (∞ → T):   f(v) = T + s·R·(v−C)/|v−C|²        [is_similarity=false]
// Project a Mat5 onto O(4,1) using a damped Schulz iteration.
//
// The standard (undamped) Schulz step  M ← M·(3I − η MᵀηM)/2  converges
// quadratically when ‖MᵀηM − η‖ < 1 but diverges for larger errors.
//
// The damped variant  M ← M·(I + α(I − η MᵀηM))  is gradient descent on
// the constraint E = MᵀηM − η and converges for ‖E‖ < 1/α.  Using an
// adaptive step α = min(0.5, 0.4/‖E‖) keeps the effective step small enough
// for any starting error while still converging quickly once close.
static Mat5 eta_orthonormalize(const Mat5 &M_in)
{
  Mat5 eta = Mat5::Identity();
  eta(4,4) = -1.0;
  Mat5 M = M_in;
  for (int iter = 0; iter < 200; iter++)
  {
    Mat5 E = M.transpose() * eta * M - eta;
    double err = E.norm();
    if (err < 1e-12) break;
    double alpha = std::min(0.5, 0.4 / err); // keeps ‖α·E‖ ≤ 0.4 < 1
    // Correct O(4,1) Schulz step: M ← M·(I − α·η·E).
    // Linearised error map: E' ≈ (1 − 2α)E, so α ∈ (0, ½] is unconditionally convergent.
    M = M * (Mat5::Identity() - alpha * eta * E);
  }
  return M;
}

//
// Extraction:
//   n_∞ = (0,0,0,−1,1)ᵀ  is the null vector for the point at infinity.
//   T  = image of ∞ = decode( M·n_∞ )
//   C  = preimage of ∞ = decode( M⁻¹·n_∞ ) = decode( η·Mᵀ·η·n_∞ )
//   sR = columns from  f(C + eᵢ) − T  (inversion) or  f(eᵢ) − T  (similarity)
static void decomposeMobius(const Mat5 &M,
                             Landscape::Set::Ball::Mobius &out)
{
  // Caller must supply an M that is already in (or close to) O(4,1).
  // See eta_orthonormalize().

  Mat5 eta = Mat5::Identity();
  eta(4,4) = -1.0;

  // Embed a 3-D point in the conformal model and transform it with M.
  auto xformPt = [&](const Eigen::Vector3d &p) -> Eigen::Vector3d {
    double r2 = p.squaredNorm();
    Vec5 P; P << p.x(), p.y(), p.z(), (1.0-r2)/2.0, (1.0+r2)/2.0;
    Vec5 Y = M * P;
    double w = Y(3) + Y(4);
    return (std::abs(w) > 1e-15) ? (Y.head<3>() / w).eval()
                                 : Eigen::Vector3d(1e15, 0.0, 0.0);
  };

  // Test whether M maps ∞ to ∞ (similarity case).
  Vec5 n_inf; n_inf << 0, 0, 0, -1, 1;
  Vec5 Y_inf = M * n_inf;
  double w_inf = Y_inf(3) + Y_inf(4);

  out.is_similarity = (std::abs(w_inf) < 1e-8);

  if (out.is_similarity)
  {
    // Similarity  f(v) = T + s·R·v.
    // T = f(0).  s·R·eᵢ = f(eᵢ) − T  for each unit basis vector.
    out.C = Eigen::Vector3d::Zero(); // unused
    out.T = xformPt(Eigen::Vector3d::Zero());
    Eigen::Matrix3d sR;
    for (int i = 0; i < 3; i++)
      sR.col(i) = xformPt(Eigen::Vector3d::Unit(i)) - out.T;
    // ‖sR‖_F = s·√3 for orthogonal R.  Average all three column norms for robustness.
    out.s = (sR.col(0).norm() + sR.col(1).norm() + sR.col(2).norm()) / 3.0;
    out.R = (out.s > 1e-15) ? Eigen::Matrix3d(sR / out.s) : Eigen::Matrix3d::Identity();
  }
  else
  {
    // Inversion  f(v) = T + s·R·(v−C)/|v−C|².
    // T = decode(M·n_∞).  C = decode(η·Mᵀ·η·n_∞).
    out.T = Y_inf.head<3>() / w_inf;

    Vec5 Z = eta * (M.transpose() * (eta * n_inf));
    double w_C = Z(3) + Z(4);
    out.C = (std::abs(w_C) > 1e-15) ? (Z.head<3>() / w_C).eval()
                                     : Eigen::Vector3d::Zero();

    // At v = C + eᵢ:  |v−C|² = 1, so f(C+eᵢ) − T = s·R·eᵢ.
    // Three test points give all three columns of s·R.
    Eigen::Matrix3d sR;
    for (int i = 0; i < 3; i++)
      sR.col(i) = xformPt(out.C + Eigen::Vector3d::Unit(i)) - out.T;
    out.s = (sR.col(0).norm() + sR.col(1).norm() + sR.col(2).norm()) / 3.0;
    out.R = (out.s > 1e-15) ? Eigen::Matrix3d(sR / out.s) : Eigen::Matrix3d::Identity();
  }
}

static bool computeMobiusTransform(Landscape::Set::Ball &src,
                                   Landscape::Set::Ball &dst)
{
  src.mobius.M = Mat5::Identity();

  int m = (int)src.type_to_set.size(); // index 0 = ball itself
  const std::string tag = "[mobius " + src.parent_set->name
                        + " ball " + std::to_string(src.type_to_set[0]) + "]";

  if (m < 2)
  {
    std::cout << tag << " identity (isolated)\n";
    return true;
  }

  Landscape::Set *src_set = src.parent_set;
  Landscape::Set *dst_set = dst.parent_set;

  Mat5 eta = Mat5::Identity();
  eta(4,4) = -1.0;

  // Build 5×m matrices A (source) and C (dest), columns = unit σ̂.
  Eigen::MatrixXd A(5, m), C(5, m);
  for (int i = 0; i < m; i++)
  {
    const auto &sb = src_set->balls[src.type_to_set[i]];
    const auto &db = dst_set->balls[dst.type_to_set[i]];
    if (sb.curvature == 0.0 || db.curvature == 0.0)
    {
      std::cerr << tag << " SKIP: plane not supported\n";
      return false;
    }
    double rs = 1.0/sb.curvature, rd = 1.0/db.curvature;
    A.col(i) = conformal_sphere(sb.dir*(sb.dist+rs), rs);
    C.col(i) = conformal_sphere(db.dir*(db.dist+rd), rd);
  }

  // Gram matrix compatibility check.  Print a warning for mismatches but
  // continue — a small mismatch (solver imprecision) gives an approximate M.
  double gram_err = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      gram_err = std::max(gram_err, std::abs(
          mink_dot(A.col(i), A.col(j)) - mink_dot(C.col(i), C.col(j))));
  if (gram_err > 1e-3)
  {
    std::cout << tag << " [Gram=" << gram_err << "] ";
    if (gram_err > 0.1)
    {
      // Print the full pair-by-pair inversive distances (= -mink_dot) for diagnosis.
      // δ(i,j) = -mink_dot(σ̂_i, σ̂_j); equals cos(intersection_angle) for tangent spheres.
      std::cout << "\n" << tag << " src/dst inversive distances (set-ball indices):\n";
      for (int i = 0; i < m; i++)
        for (int j = 0; j < i; j++)
        {
          double ds = -mink_dot(A.col(i), A.col(j));
          double dd = -mink_dot(C.col(i), C.col(j));
          int si = src.type_to_set[i], sj = src.type_to_set[j];
          int di = dst.type_to_set[i], dj = dst.type_to_set[j];
          std::cout << "    (" << si << "," << sj << ")→(" << di << "," << dj << "):"
                    << "  src δ=" << ds << "  dst δ=" << dd
                    << (std::abs(ds-dd) > 0.1 ? "  *** MISMATCH" : "") << "\n";
        }
    }
  }

  // Rank determination and M computation via thin SVD of A.
  //
  // SVD is used instead of QR because column-pivoting QR can over-estimate the
  // rank for geometrically symmetric configurations.  Example: a 4-fold
  // equatorial ring of spheres satisfies σ̂₁ + σ̂₃ = σ̂₂ + σ̂₄ exactly
  // (identical last two conformal coordinates), making A truly rank 4 even
  // though m = 5.  QR misses this; (A·Aᵀ)⁻¹ then blows up to NaN.
  //
  // With SVD the rank-deficient case is correctly routed to the η-complement
  // branch, where M = C5 · A5⁻¹ gives the identity for any self-map.
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_A(A, Eigen::ComputeFullU | Eigen::ComputeThinV);
  // Use a relative threshold of 1% of the largest singular value.
  // This correctly classifies near-zero singular values that arise from
  // exact linear dependencies in conformal space (e.g. the 4-fold equatorial
  // ring satisfies σ̂₁+σ̂₃ = σ̂₂+σ̂₄, making a 5-ball neighbourhood rank-4
  // despite having 5 columns). A tight absolute threshold (1e-8) lets these
  // near-zero values through, causing 1/sv amplification of small Gram errors.
  svd_A.setThreshold(0.01);
  int r = svd_A.rank();

  Mat5 M;
  if (r >= 5)
  {
    // Full rank (m ≥ 5): M = C · A⁺  where  A⁺ = V · S⁻¹ · Uᵀ.
    // For m = 5 this equals C · A⁻¹.  For m > 5 it is the minimum-residual
    // least-squares solution over all m columns.
    Eigen::VectorXd sinv = svd_A.singularValues().array().inverse(); // all non-zero
    M = (C * svd_A.matrixV() * sinv.asDiagonal() * svd_A.matrixU().transpose()).eval();
  }
  else
  {
    // Underdetermined (r < 5), including cases where m ≥ 5 but columns of A
    // are linearly dependent in conformal space.
    //
    // Constrained directions (k = 0..r−1): source basis = SVD left-singular
    // vectors uₖ; image = C · vₖ / sₖ.
    // Free directions (k = r..4): identity on the η-orthogonal complement
    // of span(A), which guarantees M ∈ O(4,1) on the unconstrained subspace.
    //
    // For any self-map (C = A): C·vₖ/sₖ = A·vₖ/sₖ = uₖ, so C5 = A5 → M = I.
    Eigen::MatrixXd A5(5,5), C5(5,5);
    for (int k = 0; k < r; k++)
    {
      A5.col(k) = svd_A.matrixU().col(k);
      C5.col(k) = C * svd_A.matrixV().col(k) / svd_A.singularValues()(k);
    }
    // η-complement of span(A): null space of Aᵀη.
    // For a cross-set map (src ≠ dst), C may have a *different* η-null
    // direction than A if the Gram matrices don't match exactly.  Using A's
    // null vector for both A5 and C5 only works for self-maps (C = A).
    // For cross-maps we compute C's η-null independently and let M map
    // A's null direction to C's null direction — that is, A5[:,k] = null_A
    // and C5[:,k] = null_C.  This makes Gram_η(C5) = Gram_η(A5) exactly
    // (assuming the 4D sub-Gram matrices match, which the solver enforces),
    // and therefore M = C5·A5⁻¹ ∈ O(4,1) to machine precision.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_c(A.transpose() * eta,
                                             Eigen::ComputeFullV);
    int nc = 5 - r;
    Eigen::MatrixXd null_A = svd_c.matrixV().rightCols(nc); // η-null of A
    A5.rightCols(nc) = null_A;

    // Compute the η-null space of C independently.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_cC(C.transpose() * eta,
                                              Eigen::ComputeFullV);
    Eigen::MatrixXd null_C = svd_cC.matrixV().rightCols(nc); // η-null of C
    // Align sign so that M maps null_A to null_C (not to -null_C), keeping
    // det(M) > 0 where possible.
    for (int k = 0; k < nc; k++)
      if (null_A.col(k).dot(null_C.col(k)) < 0) null_C.col(k) *= -1.0;
    C5.rightCols(nc) = null_C;

    M = C5 * A5.inverse();
  }

  // Lorentzian polarization: iterate M ← M·(3I − η MᵀηM)/2 to re-project onto
  // O(4,1) when a small Gram mismatch or numerical error has pushed M off.
  //
  // Fixed-point check: M ∈ O(4,1) ⟹ MᵀηM = η ⟹ η·Mᵀ·η·M = η² = I
  //   ⟹  M·(3I − I)/2 = M·I = M  ✓
  //
  // This is the O(4,1) analogue of the Schulz iteration for O(n): it removes
  // the η-symmetric part of the perturbation E (when M = Q(I+E), Q ∈ O(4,1))
  // with quadratic convergence, leaving only the η-skew-symmetric part.
  // (The wrong formula (3I − MᵀηM)/2 has η as a fixed point, not I.)
  //
  // Only attempt this if M is already close to O(4,1).
  double pre_err = (M.transpose() * eta * M - eta).norm();
  if (pre_err < 1.0)
  {
    for (int iter = 0; iter < 20; iter++)
    {
      Mat5 dev = M.transpose() * eta * M - eta;
      if (dev.norm() < 1e-10) break;  // generous threshold: a few ε_machine above zero
      M = M * (3.0 * Mat5::Identity() - eta * M.transpose() * eta * M) / 2.0;
    }
  }

  // Diagnostics: O(4,1) condition and full residual over all m columns.
  double oo1_err  = (M.transpose() * eta * M - eta).norm();
  double residual = (M * A - C).norm() / std::sqrt((double)m);

  std::cout << tag << "  r=" << r
            << "  O(4,1)_err=" << oo1_err
            << "  residual=" << residual;

  const double tol_ok   = 1e-4;
  const double tol_fail = 5.0;//0.1;

  // Project M onto O(4,1) via damped Schulz iteration so that the stored
  // matrix is a valid Möbius transform even when the linear solve left it
  // with residual O(4,1) error.  The projected matrix Mort is used for both
  // transformPoint and the T,C,s,R decomposition so the two are consistent.
  if (M.allFinite()) {
    const Mat5 Mort = eta_orthonormalize(M);
    double oo1_after = (Mort.transpose() * eta * Mort - eta).norm();
    std::cout << "  oo1_after=" << oo1_after;

    src.mobius.M = Mort;
    decomposeMobius(Mort, src.mobius);

    // Verify internal consistency: transformPoint (uses Mort) must match
    // transformDecomposed (uses T,C,s,R extracted from Mort).
    static const Eigen::Vector3d test_pts[] = {
      {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {0.5,0.3,-0.7}, {2,-1,3}
    };
    double decomp_err = 0.0;
    for (const auto &pt : test_pts)
    {
      auto ref  = src.mobius.transformPoint(pt);
      auto fast = src.mobius.transformDecomposed(pt);
      if (ref.allFinite() && fast.allFinite())
        decomp_err = std::max(decomp_err, (ref - fast).norm());
    }
    std::cout << "  decomp_err=" << decomp_err;
  }

  // Treat NaN as FAIL: ieee nan comparisons always return false, so guard explicitly.
  bool bad = !std::isfinite(oo1_err) || !std::isfinite(residual)
             || oo1_err > tol_fail || residual > tol_fail;
  if (bad)
  {
    std::cout << "  FAIL\n";
    return false;
  }

  std::cout << ((oo1_err > tol_ok || residual > tol_ok) ? "  APPROX\n" : "  OK\n");
  return true;
}

// ── Joint global Gauss-Seidel ─────────────────────────────────────────────────
//
// Runs a single Gauss-Seidel iteration pool that simultaneously enforces:
//   (A) Each set's own pairwise connectivity constraints (angle / distance / gap),
//       identical to the per-set applyConnectivity().
//   (B) For every ball B in set S that has a dest_set link to ball D in set T:
//       the inversive distance between every pair of type-neighbours in S must
//       equal the inversive distance between the corresponding pair in T.
//
// Inversive distance between spheres (Ci,ri) and (Cj,rj):
//   δ(i,j) = (|Ci−Cj|² − ri² − rj²) / (2 ri rj)
// This is the sole Möbius-invariant scalar for a sphere pair, so enforcing
//   δS(i,j) = δT(i,j)  for all neighbour pairs
// is necessary and sufficient for a Möbius transform to exist mapping the
// source neighbourhood to the dest neighbourhood.
//
// For each such Gram constraint the gradient of the error
//   e = δS(i,j) − δT(i,j)
// is computed w.r.t. the 6 degrees of freedom of each sphere
// (dir [2-dof tangent], dist, curvature) and a minimum-norm step is applied.
// The gradient is split evenly so that both sets are pulled toward each other.
//
void Landscape::applyConnectivity(int iterations)
{
  // ── Step 1: resolve dest_ball pointers by name ──────────────────────────
  // (mirrors the logic in matchUpDestinationBalls so we don't need it called first)
  for (auto &set : sets)
  {
    for (auto &ball : set.balls)
    {
      if (ball.dest_set == "")
      {
        ball.dest_ball = &ball; // self-map
        continue;
      }
      // Find the set with the right name.
      Set *dest_set_ptr = nullptr;
      for (auto &s : sets)
        if (s.name == ball.dest_set) { dest_set_ptr = &s; break; }
      if (!dest_set_ptr)
      {
        std::cerr << "[applyConnectivity] dest_set '" << ball.dest_set
                  << "' not found for ball in set '" << set.name << "'\n";
        continue;
      }
      if (ball.dest_ball_id >= 0 && ball.dest_ball_id < (int)dest_set_ptr->balls.size())
        ball.dest_ball = &dest_set_ptr->balls[ball.dest_ball_id];
      else
        ball.dest_ball = &dest_set_ptr->balls[0]; // fallback
    }
  }

  // ── Step 2: build all constraint pairs ──────────────────────────────────
  //
  // Each entry is one scalar constraint.  We use a tagged union approach:
  // either a within-set pair (conn order), or a cross-set Gram-matching pair.
  const double damping = 1e-10;
  const double k = 0.0; // minimum extra gap for separation constraints

  // Within-set pairs: {set index, ball i, ball j}
  struct IntraPair { int si, i, j; };
  std::vector<IntraPair> intra_pairs;
  for (int si = 0; si < (int)sets.size(); si++)
  {
    int n = (int)sets[si].balls.size();
    for (int i = 0; i < n; i++)
      for (int j = 0; j < i; j++)
        intra_pairs.push_back({si, i, j});
  }

  // Cross-set Gram pairs: for each ball B with a non-self dest_ball D,
  // emit one constraint per pair of type-neighbours of B.
  // type_to_set[ti] = set ball index; index 0 is the ball itself, ti > 0 are
  // its canonical neighbours.  addSetToTypes must be called first (main.cpp
  // already guarantees this ordering).
  //
  // We emit pairs (ti, tj) with ti < tj.  We store {ball_in_S, ball_in_T, ti, tj}.
  struct GramPair { Set::Ball *src; Set::Ball *dst; int ti, tj; };
  std::vector<GramPair> gram_pairs;
  for (auto &set : sets)
  {
    for (auto &ball : set.balls)
    {
      if (ball.dest_ball == nullptr || ball.dest_ball == &ball) continue;

      const std::vector<int> &src_nbrs = ball.type_to_set;            // ti → src set-ball idx
      const std::vector<int> &dst_nbrs = ball.dest_ball->type_to_set; // ti → dst set-ball idx

      int m = (int)std::min(src_nbrs.size(), dst_nbrs.size());
      for (int ti = 0; ti < m; ti++)
        for (int tj = 0; tj < ti; tj++)
          gram_pairs.push_back({&ball, ball.dest_ball, ti, tj});
    }
  }

  std::mt19937 rng(42);

  // ── Step 3: iterate ──────────────────────────────────────────────────────
  for (int it = 0; it < iterations; it++)
  {
    // Shuffle both pools independently to avoid ordering bias.
    std::shuffle(intra_pairs.begin(), intra_pairs.end(), rng);
    std::shuffle(gram_pairs.begin(), gram_pairs.end(), rng);

    // ── (A) Within-set constraints ─────────────────────────────────────────
    for (auto [si, i, j] : intra_pairs)
    {
      Set &set = sets[si];
      int order = set.conn(i, j);
      Set::Ball &bi = set.balls[i];
      Set::Ball &bj = set.balls[j];

      Eigen::Vector3d g_dir_i = Eigen::Vector3d::Zero(), g_dir_j = Eigen::Vector3d::Zero();
      double g_dist_i = 0, g_curv_i = 0, g_dist_j = 0, g_curv_j = 0;
      double error = 0;

      if (bi.curvature != 0.0 && bj.curvature != 0.0)
      {
        double ri = 1.0 / bi.curvature, rj = 1.0 / bj.curvature;
        Eigen::Vector3d Ci = bi.dir * (bi.dist + ri);
        Eigen::Vector3d Cj = bj.dir * (bj.dist + rj);
        Eigen::Vector3d Delta = Ci - Cj;
        double d = Delta.norm();
        if (d < 1e-12) continue;

        if (order <= 0)
        {
          double targ_d = ri + rj + (order == 0 ? k : 0.0);
          error = d - targ_d;
          if (order == 0 && error >= 0.0) continue;
          Eigen::Vector3d dddCi =  Delta / d;
          Eigen::Vector3d dddCj = -Delta / d;
          g_dist_i = dddCi.dot(bi.dir);
          g_dist_j = dddCj.dot(bj.dir);
          g_curv_i = (dddCi.dot(bi.dir) - 1.0) * (-1.0 / (bi.curvature * bi.curvature));
          g_curv_j = (dddCj.dot(bj.dir) - 1.0) * (-1.0 / (bj.curvature * bj.curvature));
          g_dir_i = (bi.dist + ri) * (dddCi - dddCi.dot(bi.dir) * bi.dir);
          g_dir_j = (bj.dist + rj) * (dddCj - dddCj.dot(bj.dir) * bj.dir);
        }
        else
        {
          double d2 = d * d;
          double cos_theta = (d2 - ri*ri - rj*rj) / (2.0 * ri * rj);
          if (cos_theta >= 1.0)
          {
            error = d - (ri + rj);
            Eigen::Vector3d dddCi =  Delta / d;
            Eigen::Vector3d dddCj = -Delta / d;
            g_dist_i = dddCi.dot(bi.dir);
            g_dist_j = dddCj.dot(bj.dir);
            g_curv_i = (dddCi.dot(bi.dir) - 1.0) * (-1.0 / (bi.curvature * bi.curvature));
            g_curv_j = (dddCj.dot(bj.dir) - 1.0) * (-1.0 / (bj.curvature * bj.curvature));
            g_dir_i = (bi.dist + ri) * (dddCi - dddCi.dot(bi.dir) * bi.dir);
            g_dir_j = (bj.dist + rj) * (dddCj - dddCj.dot(bj.dir) * bj.dir);
          }
          else if (cos_theta <= -1.0)
          {
            error = d - (ri + rj);
            Eigen::Vector3d dddCi =  Delta / d;
            Eigen::Vector3d dddCj = -Delta / d;
            g_dist_i = dddCi.dot(bi.dir);
            g_dist_j = dddCj.dot(bj.dir);
            g_curv_i = (dddCi.dot(bi.dir) - 1.0) * (-1.0 / (bi.curvature * bi.curvature));
            g_curv_j = (dddCj.dot(bj.dir) - 1.0) * (-1.0 / (bj.curvature * bj.curvature));
            g_dir_i = (bi.dist + ri) * (dddCi - dddCi.dot(bi.dir) * bi.dir);
            g_dir_j = (bj.dist + rj) * (dddCj - dddCj.dot(bj.dir) * bj.dir);
          }
          else
          {
            double theta = std::acos(cos_theta);
            double sin_theta = std::sin(theta);
            if (std::abs(sin_theta) < 1e-10) continue;
            error = pi / (double)order - theta;
            double inv_sin = 1.0 / sin_theta;
            Eigen::Vector3d dfdCi =  Delta / (ri * rj);
            Eigen::Vector3d dfdCj = -dfdCi;
            double dfdri = -(ri*ri + d2 - rj*rj) / (2.0 * ri*ri * rj);
            double dfdrj = -(rj*rj + d2 - ri*ri) / (2.0 * rj*rj * ri);
            g_dist_i = inv_sin * dfdCi.dot(bi.dir);
            g_dist_j = inv_sin * dfdCj.dot(bj.dir);
            g_curv_i = inv_sin * (dfdri + dfdCi.dot(bi.dir)) * (-1.0 / (bi.curvature * bi.curvature));
            g_curv_j = inv_sin * (dfdrj + dfdCj.dot(bj.dir)) * (-1.0 / (bj.curvature * bj.curvature));
            g_dir_i = inv_sin * (bi.dist + ri) * (dfdCi - dfdCi.dot(bi.dir) * bi.dir);
            g_dir_j = inv_sin * (bj.dist + rj) * (dfdCj - dfdCj.dot(bj.dir) * bj.dir);
          }
        }
      }
      else if (bi.curvature == 0.0 && bj.curvature == 0.0)
      {
        if (order <= 0) continue;
        double dot = bi.dir.dot(bj.dir);
        double theta = std::acos(std::clamp(dot, -1.0, 1.0));
        double sin_theta = std::sin(theta);
        if (std::abs(sin_theta) < 1e-10) continue;
        error = pi / (double)order - theta;
        double inv_sin = 1.0 / sin_theta;
        g_dir_i = inv_sin * (bj.dir - dot * bi.dir);
        g_dir_j = inv_sin * (bi.dir - dot * bj.dir);
      }
      else
      {
        const bool i_is_sphere = (bi.curvature != 0.0);
        Set::Ball &sphere = i_is_sphere ? bi : bj;
        Set::Ball &plane  = i_is_sphere ? bj : bi;
        double r = 1.0 / sphere.curvature;
        Eigen::Vector3d C = sphere.dir * (sphere.dist + r);
        double signed_dist = plane.dir.dot(C) - plane.dist;
        double cos_targ = (order <= 0) ? 1.0 : std::cos(pi / (double)order);
        double targ_dist = r * cos_targ + (order == 0 ? k : 0.0);
        error = signed_dist - targ_dist;
        if (order == 0 && error >= 0.0) continue;
        double g_dist_s = plane.dir.dot(sphere.dir);
        double g_curv_s = (plane.dir.dot(sphere.dir) - cos_targ) * (-1.0 / (sphere.curvature * sphere.curvature));
        Eigen::Vector3d g_dir_s = (sphere.dist + r) * (plane.dir - plane.dir.dot(sphere.dir) * sphere.dir);
        double g_dist_p = -1.0;
        Eigen::Vector3d g_dir_p = C - plane.dir.dot(C) * plane.dir;
        if (i_is_sphere)
        {
          g_dist_i = g_dist_s; g_curv_i = g_curv_s; g_dir_i = g_dir_s;
          g_dist_j = g_dist_p; g_curv_j = 0.0;      g_dir_j = g_dir_p;
        }
        else
        {
          g_dist_j = g_dist_s; g_curv_j = g_curv_s; g_dir_j = g_dir_s;
          g_dist_i = g_dist_p; g_curv_i = 0.0;      g_dir_i = g_dir_p;
        }
      }

      double g2 = g_dir_i.squaredNorm() + g_dist_i*g_dist_i + g_curv_i*g_curv_i
                + g_dir_j.squaredNorm() + g_dist_j*g_dist_j + g_curv_j*g_curv_j;
      double step = -error / (g2 + damping);
      bi.dir       += step * g_dir_i;  bi.dir.normalize();
      bi.dist      += step * g_dist_i;
      bi.curvature  = std::max(1e-6, bi.curvature + step * g_curv_i);
      bj.dir       += step * g_dir_j;  bj.dir.normalize();
      bj.dist      += step * g_dist_j;
      bj.curvature  = std::max(1e-6, bj.curvature + step * g_curv_j);
    }

    // ── (B) Cross-set Gram (inversive distance) constraints ────────────────
    //
    // For each pair (ti, tj) in the neighbourhood, enforce:
    //   δS(i,j) = δT(i,j)
    //   where δ(a,b) = (|Ca−Cb|² − ra² − rb²) / (2 ra rb)
    //
    // error = δS − δT
    // Gradient of δ(a,b) w.r.t. Ca: ΔCab / (ra rb)
    // Gradient of δ(a,b) w.r.t. ra: −(ra² + |ΔC|² − rb²) / (2 ra² rb)
    //                              = −(δ(a,b)/ra + 1/rb)   [alternative form]
    //
    // Only the source side is updated; the dest geometry is fixed by its own
    // intra constraints and updating it would cause oscillation.
    for (auto &gp : gram_pairs)
    {
      Set::Ball &src_ball = *gp.src;
      Set::Ball &dst_ball = *gp.dst;
      const auto &src_n2s = src_ball.type_to_set; // may be empty if addSetToTypes not called
      const auto &dst_n2s = dst_ball.type_to_set;

      // Resolve set-ball indices for this (ti, tj) pair.
      int si_idx = src_n2s[gp.ti], sj_idx = src_n2s[gp.tj];
      int di_idx = dst_n2s[gp.ti], dj_idx = dst_n2s[gp.tj];

      Set &src_set = *src_ball.parent_set;
      Set &dst_set = *dst_ball.parent_set;
      Set::Ball &sA = src_set.balls[si_idx];
      Set::Ball &sB = src_set.balls[sj_idx];
      Set::Ball &dA = dst_set.balls[di_idx];
      Set::Ball &dB = dst_set.balls[dj_idx];

      // Only handle sphere-sphere pairs for now (planes have infinite inversive dist).
      if (sA.curvature == 0.0 || sB.curvature == 0.0 ||
          dA.curvature == 0.0 || dB.curvature == 0.0) continue;

      double rSA = 1.0/sA.curvature, rSB = 1.0/sB.curvature;
      double rDA = 1.0/dA.curvature, rDB = 1.0/dB.curvature;
      Eigen::Vector3d CSA = sA.dir*(sA.dist+rSA), CSB = sB.dir*(sB.dist+rSB);
      Eigen::Vector3d CDA = dA.dir*(dA.dist+rDA), CDB = dB.dir*(dB.dist+rDB);

      auto inv_dist = [](const Eigen::Vector3d &Ca, double ra,
                         const Eigen::Vector3d &Cb, double rb)
      {
        return ((Ca-Cb).squaredNorm() - ra*ra - rb*rb) / (2.0*ra*rb);
      };
      double delta_S = inv_dist(CSA, rSA, CSB, rSB);
      double delta_T = inv_dist(CDA, rDA, CDB, rDB);
      double error   = delta_S - delta_T;

      // Gradients of δS w.r.t. the four source spheres:
      //   dδ/dCa = (Ca−Cb)/(ra rb)
      //   dδ/dra = −(ra² + |ΔC|² − rb²) / (2 ra² rb)  =  −(δ+1/rb) / ra
      //            but computed directly to avoid cancellation.
      double d2S = (CSA-CSB).squaredNorm();
      double inv_rSArSB = 1.0/(rSA*rSB);
      Eigen::Vector3d gCdCSA =  (CSA-CSB) * inv_rSArSB;
      Eigen::Vector3d gCdCSB = -(CSA-CSB) * inv_rSArSB;
      double gCdrSA = -(rSA*rSA + d2S - rSB*rSB) / (2.0*rSA*rSA*rSB);
      double gCdrSB = -(rSB*rSB + d2S - rSA*rSA) / (2.0*rSB*rSB*rSA);

      // Chain to (dir, dist, curvature): dC/d(dist) = dir, dC/d(dir_tangent) = (dist+r)*I_tang,
      // dr/d(curvature) = -1/curvature².
      auto grad_for = [](const Set::Ball &b, double r,
                         const Eigen::Vector3d &gC, double gr)
        -> std::tuple<Eigen::Vector3d, double, double>
      {
        Eigen::Vector3d g_dir  = (b.dist+r)*(gC - gC.dot(b.dir)*b.dir);
        double          g_dist = gC.dot(b.dir);
        double          g_curv = (gr + gC.dot(b.dir)) * (-1.0/(b.curvature*b.curvature));
        return {g_dir, g_dist, g_curv};
      };

      auto [gDirSA, gDistSA, gCurvSA] = grad_for(sA, rSA, gCdCSA, gCdrSA);
      auto [gDirSB, gDistSB, gCurvSB] = grad_for(sB, rSB, gCdCSB, gCdrSB);

      // Gradients of δT are the negatives of δS gradients (error = δS − δT),
      // but applied to the dest spheres dA, dB.
      double d2T = (CDA-CDB).squaredNorm();
      double inv_rDArDB = 1.0/(rDA*rDB);
      Eigen::Vector3d gCdCDA =  (CDA-CDB) * inv_rDArDB;
      Eigen::Vector3d gCdCDB = -(CDA-CDB) * inv_rDArDB;
      double gCdrDA = -(rDA*rDA + d2T - rDB*rDB) / (2.0*rDA*rDA*rDB);
      double gCdrDB = -(rDB*rDB + d2T - rDA*rDA) / (2.0*rDB*rDB*rDA);

      auto [gDirDA, gDistDA, gCurvDA] = grad_for(dA, rDA, gCdCDA, gCdrDA);
      auto [gDirDB, gDistDB, gCurvDB] = grad_for(dB, rDB, gCdCDB, gCdrDB);

      double g2 = gDirSA.squaredNorm() + gDistSA*gDistSA + gCurvSA*gCurvSA
                + gDirSB.squaredNorm() + gDistSB*gDistSB + gCurvSB*gCurvSB;
                // Note: only the *source* side drives the step.  The dest set's
                // geometry is fully determined by its own intra constraints; if we
                // also pushed dest balls here they would fight back (oscillation).
      double step = -error / (g2 + damping);

      // Source side only: pull shell-shell (src) toward the dest Gram matrix.
      sA.dir += step * gDirSA;  sA.dir.normalize();
      sA.dist += step*gDistSA;
      sA.curvature = std::max(1e-6, sA.curvature + step*gCurvSA);
      sB.dir += step * gDirSB;  sB.dir.normalize();
      sB.dist += step*gDistSB;
      sB.curvature = std::max(1e-6, sB.curvature + step*gCurvSB);
    }
  }
}

void Landscape::matchUpDestinationBalls()
{
  for (auto &set: sets)
  {
    for (auto &ball: set.balls)
    {
      if (ball.dest_set == "") // standard recursion
      {
        ball.dest_ball = &ball;
        computeMobiusTransform(ball, ball);
        continue;
      }
      const Type *type = ball.type;
      for (auto &type_ball: type->balls) // for every ball of this type
      {
        if (type_ball->parent_set->name == ball.dest_set) // check it matches the set name
        {
          Set *par = type_ball->parent_set;
          if (ball.dest_ball_id != -1) // if we have a specific id within the set
          {
            if (&par->balls[ball.dest_ball_id] == type_ball) // then check that id is this ball
            {
              ball.dest_ball = type_ball;
              break;
            }
          }
          else // otherwise just use the first ball that matches the set name
          {
            ball.dest_ball = type_ball;
            break;
          }
        }
      }
      if (ball.dest_ball != nullptr)
        computeMobiusTransform(ball, *ball.dest_ball);
      else
        std::cerr << "[matchUp] no dest_ball found for ball in set '" << set.name << "'\n";
    }
  }
}