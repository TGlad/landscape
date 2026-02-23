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
    std::cout << tag << " [Gram=" << gram_err << "] ";

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
  svd_A.setThreshold(1e-8);
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
    // η-complement = null space of Aᵀ η (right null space of the m×5 matrix).
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_c(A.transpose() * eta,
                                             Eigen::ComputeFullV);
    int nc = 5 - r;
    Eigen::MatrixXd null_vecs = svd_c.matrixV().rightCols(nc); // 5×nc
    A5.rightCols(nc) = null_vecs;
    C5.rightCols(nc) = null_vecs; // identity on the free complement
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
  const double tol_fail = 0.1;

  // Treat NaN as FAIL: ieee nan comparisons always return false, so guard explicitly.
  bool bad = !std::isfinite(oo1_err) || !std::isfinite(residual)
             || oo1_err > tol_fail || residual > tol_fail;
  if (bad)
  {
    std::cout << "  FAIL\n";
    return false;
  }

  std::cout << ((oo1_err > tol_ok || residual > tol_ok) ? "  APPROX\n" : "  OK\n");
  src.mobius.M = M;
  return true;
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