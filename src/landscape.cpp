#include "landscape.h"
#include <numeric>
#include <random>

static const double pi = std::acos(-1.0);

void Landscape::Set::applyConnectivity()
{
  // the state is 
  // Eigen::Vector3d dir;
  //  double dist;
  //  double curvature;
  // on n balls   
  // = 5 * n unknowns
  // there are up to n(n-1)/2 constraints

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

void Landscape::addSetToTypes(Set &set)
{
  for (int i = 0; i<(int)set.balls.size(); i++)
  { 
    Type new_type(set.conn, i);
    bool found = false;
    for (int j = 0; j<types.size(); j++)
    {
      if (types[j].conn == new_type.conn)
      {
        // new_type already exists, so add it in.
        types[j].balls.push_back(&set.balls[i]);
        set.balls[i].type = &types[j];
        found = true;
        break;
      }
    }
    if (!found)
    {
      types.push_back(new_type);
      types.back().balls.push_back(&set.balls[i]);
      set.balls[i].type = &types.back();
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
        continue;
      }
      const Type *type = ball.type;
      Set::Ball *dest_ball = nullptr;
      bool found_ball = false;
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
    }
  }
}