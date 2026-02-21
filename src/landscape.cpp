#include "landscape.h"

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

  // gradient descent with SOR, applies each constraint repeatedly
  for (int it = 0; it<10; it++)
  {
    for (int i = 0; i<(int)balls.size(); i++)
    {
      for (int j = 0; j<i; j++)
      {
        if (conn(i,j) == 0) // actually we need to ensure it stays out, which is a unilateral constraint
          continue;
        int order = conn(i,j);
        double targ_angle = order==-1 ? 0.0 : pi/(double)order;
        double angle = [&]() -> double 
        {
          const Ball &bi = balls[i];
          const Ball &bj = balls[j];
          if (bi.curvature != 0.0 && bj.curvature != 0.0) 
          {
            // sphere-sphere: cos θ = (d² - ri² - rj²) / (2 ri rj)
            double ri = 1.0 / bi.curvature, rj = 1.0 / bj.curvature;
            Eigen::Vector3d Ci = bi.dir * (bi.dist + ri);
            Eigen::Vector3d Cj = bj.dir * (bj.dist + rj);
            double d2 = (Ci - Cj).squaredNorm();
            double cos_theta = (d2 - ri*ri - rj*rj) / (2.0*ri*rj);
            if (cos_theta > 1.0) 
              return std::numeric_limits<double>::quiet_NaN(); // non-intersecting (too far apart)
            if (cos_theta < -1.0) 
              return std::numeric_limits<double>::quiet_NaN(); // non-intersecting (one inside the other)
            return std::acos(cos_theta);
          } 
          else if (bi.curvature == 0.0 && bj.curvature == 0.0) 
          {
            // plane-plane: cos θ = ni · nj
            return std::acos(std::clamp(bi.dir.dot(bj.dir), -1.0, 1.0));
          } 
          else 
          {
            // sphere-plane: cos θ = (n_plane · C_sphere - dist_plane) / r_sphere
            const Ball &sphere = (bi.curvature != 0.0) ? bi : bj;
            const Ball &plane  = (bi.curvature != 0.0) ? bj : bi;
            double r = 1.0 / sphere.curvature;
            Eigen::Vector3d C = sphere.dir * (sphere.dist + r);
            return std::acos(std::clamp((plane.dir.dot(C) - plane.dist) / r, -1.0, 1.0));
          }
        }();
        if (std::isnan(angle))
          continue; // spheres don't intersect, no constraint to apply
        double error = targ_angle - angle;

        /*
        double dAngle_dDist = ;
        Eigen;:Vector3d dAngle_dDir = ;
        double dAngle_dCurv = ;*/
      }
    }
  }
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