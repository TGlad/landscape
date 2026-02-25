#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

void Landscape::Set::Ball::initSphere(const Eigen::Vector3d &p, double rad, bool fix)
{
  dir = p.normalized();
  dist = p.norm() - rad;
  curvature = 1.0/rad;
  is_fixed = fix;
}
void Landscape::Set::Ball::initPlane(const Eigen::Vector3d &normal, double d, bool fix)
{
  dir = normal.normalized();
  dist = d;
  curvature = 0.0;
  is_fixed = fix;
}

Eigen::Vector3d Landscape::Set::Ball::Mobius::transformPoint(const Eigen::Vector3d &p) const
{
  double r2 = p.squaredNorm();
  Eigen::Matrix<double,5,1> P;
  P << p.x(), p.y(), p.z(), (1.0-r2)/2.0, (1.0+r2)/2.0;
  auto Y = M * P;
  double w = Y(3) + Y(4);
  return (std::abs(w) > 1e-15) ? Eigen::Vector3d(Y(0)/w, Y(1)/w, Y(2)/w): Eigen::Vector3d::Zero();
}

// Apply the T,C,s,R decomposition (GLSL path) to a point.
//   Similarity (is_similarity=true):  f(v) = T + s*R*v
//   Inversion  (is_similarity=false): f(v) = T + s*R*(v-C)/|v-C|²
Eigen::Vector3d Landscape::Set::Ball::Mobius::transformDecomposed(const Eigen::Vector3d &p) const
{
  if (is_similarity)
    return T + s * (R * p);
  Eigen::Vector3d w = p - C;
  double d2 = w.squaredNorm();
  if (d2 < 1e-30) return Eigen::Vector3d(1e15, 0.0, 0.0);
  return T + (s / d2) * (R * w);
}

// Returns (centre, radius) of the transformed sphere.
std::pair<Eigen::Vector3d,double> Landscape::Set::Ball::Mobius::transformSphere(
    const Eigen::Vector3d &c, double r) const
{
  double c2 = c.squaredNorm();
  Eigen::Matrix<double,5,1> sv;
  sv << c.x(), c.y(), c.z(),
        (1.0-c2+r*r)/2.0, (1.0+c2-r*r)/2.0;
  auto Y = M * (sv / r);
  double inv_r = Y(3) + Y(4);
  if (std::abs(inv_r) < 1e-15) return {{0,0,0}, 0};
  return { Y.head<3>() / inv_r, 1.0/inv_r };
}