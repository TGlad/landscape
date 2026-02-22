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

// Returns false if the Möbius transform could not be computed.
// Tries a pure similarity first; if radii scale inconsistently, falls back to
// an inversion in the source ball followed by a similarity.
static bool computeMobiusTransform(Landscape::Set::Ball &src, Landscape::Set::Ball &dst)
{
  using Vec3 = Eigen::Vector3d;
  using Mat3 = Eigen::Matrix3d;
  using Sph  = std::pair<Vec3, double>; // (centre, radius)

  // Build canonical sphere list: index 0 = the ball itself, 1..m-1 = canonical neighbours.
  auto buildSpheres = [](const Landscape::Set::Ball &b) -> std::vector<Sph>
  {
    std::vector<Sph> sph;
    const Landscape::Set *ps = b.parent_set;
    for (int si : b.type_to_set)
    {
      const Landscape::Set::Ball &nb = ps->balls[si];
      if (nb.curvature == 0.0) return {}; // planes not yet supported
      double r = 1.0 / nb.curvature;
      sph.push_back({ nb.dir * (nb.dist + r), r });
    }
    return sph;
  };

  auto src_sph = buildSpheres(src);
  auto dst_sph = buildSpheres(dst);
  if (src_sph.empty() || dst_sph.empty()) return false;
  int m = (int)src_sph.size();

  const std::string tag = "[mobius " + src.parent_set->name
                        + " ball " + std::to_string(src.type_to_set[0]) + "]";

  // Kabsch SVD: find best proper rotation R s.t. R*(from_i - from_0) ≈ (to_i - to_0)/s.
  auto kabsch = [&](const std::vector<Sph> &from, const std::vector<Sph> &to, double s) -> Mat3
  {
    if (m <= 1) return Mat3::Identity();
    Mat3 H = Mat3::Zero();
    for (int i = 1; i < m; i++)
    {
      Vec3 A = from[i].first - from[0].first;
      Vec3 B = (to[i].first  - to[0].first) / s;
      H += A * B.transpose();
    }
    Eigen::JacobiSVD<Mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Mat3 U = svd.matrixU(), V = svd.matrixV();
    Mat3 D = Mat3::Identity();
    D(2,2) = (V * U.transpose()).determinant() > 0 ? 1.0 : -1.0;
    return V * D * U.transpose();
  };

  // RMS residual of mapped centres.
  auto residual = [&](const std::vector<Sph> &from, const std::vector<Sph> &to,
                      double s, const Mat3 &R, const Vec3 &t) -> double
  {
    double e2 = 0;
    for (int i = 0; i < m; i++)
      e2 += (s * R * from[i].first + t - to[i].first).squaredNorm();
    return std::sqrt(e2 / m);
  };

  // ---- Try similarity (flip = false) ----
  {
    double s = 0;
    for (int i = 0; i < m; i++) s += dst_sph[i].second / src_sph[i].second;
    s /= m;
    bool ok = true;
    for (int i = 0; i < m; i++)
      if (std::abs(dst_sph[i].second / src_sph[i].second - s) > 1e-4 * s) { ok = false; break; }
    if (ok)
    {
      Mat3 R = kabsch(src_sph, dst_sph, s);
      Vec3 t = dst_sph[0].first - s * R * src_sph[0].first;
      double err = residual(src_sph, dst_sph, s, R, t);
      std::cout << tag << " similarity  s=" << s << "  residual=" << err << "\n";
      if (err < 1e-4)
      {
        src.mobius = { Vec3::Zero(), t, R, s, false };
        return true;
      }
    }
  }

  // ---- Try inversion in src ball, then similarity (flip = true) ----
  {
    const Vec3 &C0 = src_sph[0].first;
    double rho2    = src_sph[0].second * src_sph[0].second;

    std::vector<Sph> inv(m);
    bool ok = true;
    for (int i = 0; i < m; i++)
    {
      Vec3   dv = src_sph[i].first - C0;
      double ri = src_sph[i].second;
      double D  = dv.squaredNorm() - ri * ri;
      if (std::abs(D) < 1e-10) { ok = false; break; }
      double k = rho2 / D;
      inv[i] = { C0 + k * dv, std::abs(k) * ri };
    }

    if (ok)
    {
      double s = 0;
      for (int i = 0; i < m; i++) s += dst_sph[i].second / inv[i].second;
      s /= m;
      bool scale_ok = true;
      for (int i = 0; i < m; i++)
        if (std::abs(dst_sph[i].second / inv[i].second - s) > 1e-4 * s) { scale_ok = false; break; }

      if (scale_ok)
      {
        Mat3 R = kabsch(inv, dst_sph, s);
        Vec3 t = dst_sph[0].first - s * R * inv[0].first;
        double err = residual(inv, dst_sph, s, R, t);
        std::cout << tag << " inversion+similarity  rho2=" << rho2
                  << "  s=" << s << "  residual=" << err << "\n";
        if (err < 1e-4)
        {
          // M(p) = R * (s*rho2 * (p-C0)/|p-C0|²) + t
          // struct: rotation*(scale*(p-center)/|p-center|²) + translation
          src.mobius = { C0, t, R, s * rho2, true };
          return true;
        }
      }
    }
  }

  std::cerr << tag << " FAIL: could not compute Möbius transform\n";
  return false;
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