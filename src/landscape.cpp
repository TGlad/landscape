#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

static const double pi = std::acos(-1.0);

#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <functional>
#include <limits>
#include <Eigen/Dense>

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

void Landscape::verifyConnectivity()
{
  for (auto &set: sets)
    set.verifyConnectivity();
}

void Landscape::calculateLeaves()
{
  for (auto &set: sets)
  {
    if (set.leaf_ball_set.size() > 0)
      set.calculateLeafBalls();
    for (int i = 0; i<(int)set.leaf_ball_ids.size(); i++)
    {
      auto &lbi = set.leaf_ball_ids[i];
      if (lbi[3] == -1)
        set.calculateLeafBall(lbi[0], lbi[1], lbi[2], set.leaf_ball_scales[i]);
      else 
        set.calculateLeafBall(lbi[0], lbi[1], lbi[2], lbi[3]);
    }
  }
}
void Landscape::addSetsToTypes()
{
  for (auto &set: sets)
  {
    addSetToTypes(set);
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

static Mat5 eta_orthonormalize(const Mat5 &M_in);

static double clampSignedRadius(double r)
{
  constexpr double eps = 1e-6;
  if (std::abs(r) < eps)
    return (r < 0.0) ? -eps : eps;
  return r;
}

static double mink_dot(const Vec5 &a, const Vec5 &b)
{
  return a.head<4>().dot(b.head<4>()) - a(4)*b(4);
}

static Vec5 conformal_sphere(const Eigen::Vector3d &c, double r)
{
  r = clampSignedRadius(r);
  double c2 = c.squaredNorm();
  Vec5 v;
  v << c.x(), c.y(), c.z(), (1.0-c2+r*r)/2.0, (1.0+c2-r*r)/2.0;
  return v / r; // unit spacelike: v·η·v = 1
}

static Vec5 conformal_plane(const Eigen::Vector3d &n, double d)
{
  Vec5 v;
  v << n.x(), n.y(), n.z(), -d, d;
  return v; // unit spacelike when |n|=1: v·η·v = 1
}

static Vec5 conformal_ball(const Landscape::Set::Ball &b)
{
  return conformal_sphere(b.centre, b.radius);
}

static bool fitMobiusFromPairs(const Landscape::Set &src_set,
                               const Landscape::Set &dst_set,
                               const std::vector<std::pair<int,int>> &pairs,
                               Mat5 &M_out,
                               double &residual_out)
{
  int m = (int)pairs.size();
  if (m < 2) return false;

  Eigen::MatrixXd A(5, m), C(5, m);
  for (int i = 0; i < m; i++)
  {
    int si = pairs[i].first;
    int di = pairs[i].second;
    if (si < 0 || si >= (int)src_set.balls.size() || di < 0 || di >= (int)dst_set.balls.size())
      return false;
    A.col(i) = conformal_ball(src_set.balls[si]);
    C.col(i) = conformal_ball(dst_set.balls[di]);
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> svd_A(A, Eigen::ComputeFullU | Eigen::ComputeThinV);
  svd_A.setThreshold(0.01);
  int r = svd_A.rank();
  if (r <= 0) return false;

  Eigen::VectorXd sinv = Eigen::VectorXd::Zero(svd_A.singularValues().size());
  for (int i = 0; i < sinv.size(); i++)
    if (svd_A.singularValues()(i) > 0.0)
      sinv(i) = 1.0 / svd_A.singularValues()(i);

  Mat5 M = (C * svd_A.matrixV() * sinv.asDiagonal() * svd_A.matrixU().transpose()).eval();
  if (!M.allFinite()) return false;

  M = eta_orthonormalize(M);
  if (!M.allFinite()) return false;

  residual_out = (M * A - C).norm() / std::sqrt((double)m);
  M_out = M;
  return std::isfinite(residual_out);
}

static std::vector<int> identityTypeMap(const Landscape::Set::Ball &src,
                                        const Landscape::Set::Ball &dst)
{
  int m = (int)std::min(src.type_to_set.size(), dst.type_to_set.size());
  std::vector<int> map(m);
  for (int i = 0; i < m; i++) map[i] = i;
  return map;
}

static int mapSetBallThroughTypeMap(const Landscape::Set::Ball &src,
                                    const Landscape::Set::Ball &dst,
                                    const std::vector<int> &dst_ti_for_src_ti,
                                    int src_set_ball)
{
  if (src_set_ball < 0 || src_set_ball >= (int)src.set_to_type.size()) return -1;
  int src_ti = src.set_to_type[src_set_ball];
  if (src_ti < 0 || src_ti >= (int)dst_ti_for_src_ti.size()) return -1;
  int dst_ti = dst_ti_for_src_ti[src_ti];
  if (dst_ti < 0 || dst_ti >= (int)dst.type_to_set.size()) return -1;
  return dst.type_to_set[dst_ti];
}

static std::vector<int> orderedFanByCentroid(const Landscape::Set &set, int center_ball)
{
  std::vector<int> neigh;
  if (center_ball < 0 || center_ball >= (int)set.balls.size()) return neigh;

  for (int i = 0; i < (int)set.balls.size(); i++)
    if (i != center_ball && set.conn(center_ball, i) > 0)
      neigh.push_back(i);

  if (neigh.size() <= 1) return neigh;

  // Topology-first ring ordering: walk the induced ring graph among neighbours
  // of center_ball using edges with conn > 0.
  std::unordered_map<int, std::vector<int>> ring_adj;
  for (int a : neigh) ring_adj[a] = {};
  for (int i = 0; i < (int)neigh.size(); i++)
  {
    for (int j = i + 1; j < (int)neigh.size(); j++)
    {
      int a = neigh[i], b = neigh[j];
      if (set.conn(a, b) > 0)
      {
        ring_adj[a].push_back(b);
        ring_adj[b].push_back(a);
      }
    }
  }

  bool ring_like = true;
  for (int a : neigh)
    if ((int)ring_adj[a].size() != 2)
      ring_like = false;

  std::vector<int> ordered;
  if (ring_like)
  {
    int start = *std::min_element(neigh.begin(), neigh.end());
    int prev = -1;
    int cur = start;
    for (int step = 0; step < (int)neigh.size(); step++)
    {
      ordered.push_back(cur);
      int n0 = ring_adj[cur][0];
      int n1 = ring_adj[cur][1];
      int next = (n0 == prev) ? n1 : n0;
      prev = cur;
      cur = next;
      if (cur == start && step + 1 < (int)neigh.size())
      {
        ring_like = false;
        break;
      }
    }
  }

  if (!ring_like || (int)ordered.size() != (int)neigh.size())
  {
    // Fallback for non-ring neighbourhoods: geometric angular ordering.
    ordered.clear();

    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    for (const auto &b : set.balls) centroid += b.centre;
    centroid /= (double)std::max(1, (int)set.balls.size());

    const Eigen::Vector3d C0 = set.balls[center_ball].centre;
    Eigen::Vector3d axis = C0 - centroid;
    if (axis.squaredNorm() < 1e-20) axis = Eigen::Vector3d(0, 0, 1);
    axis.normalize();

    auto tangentComponent = [&](int bi) {
      Eigen::Vector3d v = set.balls[bi].centre - C0;
      return v - v.dot(axis) * axis;
    };

    Eigen::Vector3d u = tangentComponent(neigh[0]);
    if (u.squaredNorm() < 1e-20)
    {
      for (int i = 1; i < (int)neigh.size(); i++)
      {
        u = tangentComponent(neigh[i]);
        if (u.squaredNorm() >= 1e-20) break;
      }
    }
    if (u.squaredNorm() < 1e-20)
    {
      std::sort(neigh.begin(), neigh.end());
      ordered = neigh;
    }
    else
    {
      u.normalize();
      Eigen::Vector3d v_axis = axis.cross(u);
      if (v_axis.squaredNorm() < 1e-20)
      {
        std::sort(neigh.begin(), neigh.end());
        ordered = neigh;
      }
      else
      {
        v_axis.normalize();
        std::vector<std::pair<double, int>> ang;
        ang.reserve(neigh.size());
        for (int nb : neigh)
        {
          Eigen::Vector3d t = tangentComponent(nb);
          if (t.squaredNorm() < 1e-20)
          {
            ang.push_back({0.0, nb});
            continue;
          }
          t.normalize();
          double a = std::atan2(t.dot(v_axis), t.dot(u));
          ang.push_back({a, nb});
        }

        std::sort(ang.begin(), ang.end(), [](const auto &A, const auto &B) {
          if (A.first == B.first) return A.second < B.second;
          return A.first < B.first;
        });
        ordered.reserve(ang.size());
        for (const auto &p : ang) ordered.push_back(p.second);
      }
    }
  }

  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (const auto &b : set.balls) centroid += b.centre;
  centroid /= (double)std::max(1, (int)set.balls.size());

  const Eigen::Vector3d C0 = set.balls[center_ball].centre;
  // Handedness convention: mean ring normal should align with center->centroid.
  Eigen::Vector3d outward = centroid - C0;
  Eigen::Vector3d mean_n = Eigen::Vector3d::Zero();
  for (int i = 0; i < (int)ordered.size(); i++)
  {
    int a = ordered[i];
    int b = ordered[(i + 1) % ordered.size()];
    Eigen::Vector3d va = set.balls[a].centre - C0;
    Eigen::Vector3d vb = set.balls[b].centre - C0;
    mean_n += va.cross(vb);
  }

  if (mean_n.squaredNorm() > 1e-20 && outward.squaredNorm() > 1e-20)
  {
    if (mean_n.dot(outward) < 0.0)
      std::reverse(ordered.begin(), ordered.end());
  }

  return ordered;
}

static std::vector<int> buildFanTypeMap(const Landscape::Set::Ball &src,
                                        const Landscape::Set::Ball &dst,
                                        int src_anchor_set_ball,
                                        int dst_anchor_set_ball)
{
  std::vector<int> map = identityTypeMap(src, dst);
  if (map.empty()) return map;
  map[0] = 0;

  const int src0 = src.type_to_set[0];
  const int dst0 = dst.type_to_set[0];
  auto src_fan = orderedFanByCentroid(*src.parent_set, src0);
  auto dst_fan = orderedFanByCentroid(*dst.parent_set, dst0);
  if (src_fan.empty() || src_fan.size() != dst_fan.size()) return map;

  int src_anchor_pos = -1, dst_anchor_pos = -1;
  for (int i = 0; i < (int)src_fan.size(); i++) if (src_fan[i] == src_anchor_set_ball) src_anchor_pos = i;
  for (int i = 0; i < (int)dst_fan.size(); i++) if (dst_fan[i] == dst_anchor_set_ball) dst_anchor_pos = i;
  if (src_anchor_pos < 0 || dst_anchor_pos < 0) return map;

  int N = (int)src_fan.size();
  for (int i = 0; i < N; i++)
  {
    int s_ball = src_fan[i];
    int d_ball = dst_fan[(i - src_anchor_pos + dst_anchor_pos + N) % N];
    if (s_ball < 0 || s_ball >= (int)src.set_to_type.size()) continue;
    if (d_ball < 0 || d_ball >= (int)dst.set_to_type.size()) continue;
    int s_ti = src.set_to_type[s_ball];
    int d_ti = dst.set_to_type[d_ball];
    if (s_ti < 0 || d_ti < 0) continue;
    if (s_ti < (int)map.size() && d_ti < (int)dst.type_to_set.size())
      map[s_ti] = d_ti;
  }

  return map;
}

static void printMappedFanAdjacencyDiagnostics(const Landscape::Set::Ball &src,
                                               const Landscape::Set::Ball &dst,
                                               const std::vector<int> &dst_ti_for_src_ti,
                                               const std::string &label)
{
  if (src.parent_set == nullptr || dst.parent_set == nullptr) return;
  if (src.type_to_set.empty() || dst.type_to_set.empty()) return;

  const int src0 = src.type_to_set[0];
  const int dst0 = dst.type_to_set[0];
  auto src_fan = orderedFanByCentroid(*src.parent_set, src0);
  auto dst_fan = orderedFanByCentroid(*dst.parent_set, dst0);
  if (src_fan.empty())
  {
    std::cout << "[mapped-fan] " << label << " no source fan\n";
    return;
  }

  std::vector<int> mapped_fan;
  mapped_fan.reserve(src_fan.size());
  bool all_mapped = true;
  for (int s_ball : src_fan)
  {
    int d_ball = mapSetBallThroughTypeMap(src, dst, dst_ti_for_src_ti, s_ball);
    if (d_ball < 0) all_mapped = false;
    mapped_fan.push_back(d_ball);
  }

  std::cout << "[mapped-fan] " << label
            << " src_center=" << src0
            << " dst_center=" << dst0
            << " src:";
  for (int b : src_fan) std::cout << " " << b;
  std::cout << " mapped:";
  for (int b : mapped_fan) std::cout << " " << b;
  std::cout << "\n";

  if (!all_mapped)
  {
    std::cout << "[mapped-fan] " << label << " incomplete mapping\n";
    return;
  }

  bool center_ok = true;
  bool ring_ok = true;
  bool unique_ok = true;
  std::set<int> used;
  const int N = (int)mapped_fan.size();
  for (int i = 0; i < N; i++)
  {
    int a = mapped_fan[i];
    int b = mapped_fan[(i + 1) % N];
    if (!used.insert(a).second) unique_ok = false;
    if (a < 0 || a >= (int)dst.parent_set->balls.size())
    {
      center_ok = false;
      ring_ok = false;
      continue;
    }
    if (dst.parent_set->conn(dst0, a) <= 0) center_ok = false;
    if (b < 0 || b >= (int)dst.parent_set->balls.size() || dst.parent_set->conn(a, b) <= 0)
      ring_ok = false;
  }

  auto isCyclicShift = [&](const std::vector<int> &ref, const std::vector<int> &seq) {
    if (ref.size() != seq.size()) return false;
    int n = (int)ref.size();
    for (int sh = 0; sh < n; sh++)
    {
      bool ok = true;
      for (int i = 0; i < n; i++)
      {
        if (seq[i] != ref[(i + sh) % n])
        {
          ok = false;
          break;
        }
      }
      if (ok) return true;
    }
    return false;
  };

  bool cyclic_ok = isCyclicShift(dst_fan, mapped_fan);

  std::cout << "[mapped-fan] " << label
            << " center-adjacency=" << (center_ok ? "OK" : "BAD")
            << " ring-adjacency=" << (ring_ok ? "OK" : "BAD")
            << " unique=" << (unique_ok ? "OK" : "BAD")
            << " cyclic-with-dst-fan=" << (cyclic_ok ? "OK" : "BAD")
            << "\n";
}

static std::vector<int> buildInheritedTypeMapBySetIndex(
    const Landscape::Set::Ball &src,
    const Landscape::Set::Ball &dst,
    const Landscape::Set::Ball &owner,
    const std::vector<int> &owner_map)
{
  std::vector<int> map = identityTypeMap(src, dst);
  for (int sti = 0; sti < (int)map.size(); sti++)
  {
    if (sti < 0 || sti >= (int)src.type_to_set.size()) continue;
    const int src_set_ball = src.type_to_set[sti];
    const int dst_set_ball = mapSetBallThroughTypeMap(owner, *owner.dest_ball, owner_map, src_set_ball);
    if (dst_set_ball < 0 || dst_set_ball >= (int)dst.set_to_type.size()) continue;
    const int dst_ti = dst.set_to_type[dst_set_ball];
    if (dst_ti < 0 || dst_ti >= (int)dst.type_to_set.size()) continue;
    map[sti] = dst_ti;
  }
  return map;
}

static std::vector<int> findBestNeighbourTypeMap(const Landscape::Set::Ball &src,
                                                 const Landscape::Set::Ball &dst)
{
  std::vector<int> best_map = identityTypeMap(src, dst);
  int m = (int)best_map.size();
  if (m < 2) return best_map;

  const int src0 = src.type_to_set[0];
  const int dst0 = dst.type_to_set[0];
  auto src_fan = orderedFanByCentroid(*src.parent_set, src0);
  auto dst_fan = orderedFanByCentroid(*dst.parent_set, dst0);

  auto printFanOrderingDiagnostics = [&](const Landscape::Set &set,
                                         int center_ball,
                                         const std::vector<int> &fan,
                                         const std::string &label)
  {
    std::cout << "[fan-order] " << label
              << " center=" << center_ball
              << " fan:";
    for (int b : fan) std::cout << " " << b;
    std::cout << "\n";

    if (fan.size() < 2)
    {
      std::cout << "[fan-order] " << label << " insufficient fan size\n";
      return;
    }

    bool all_ok = true;
    int N = (int)fan.size();
    for (int i = 0; i < N; i++)
    {
      int a = fan[i];
      int b = fan[(i + 1) % N];
      int c = set.conn(a, b);
      bool ok = c > 0;
      all_ok = all_ok && ok;
      std::cout << "[fan-order] " << label
                << " edge " << a << "-" << b
                << " conn=" << c
                << (ok ? " OK" : " BAD")
                << "\n";
    }
    std::cout << "[fan-order] " << label
              << " ring-adjacency=" << (all_ok ? "OK" : "BAD")
              << "\n";
  };
  if (src_fan.empty() || src_fan.size() != dst_fan.size())
    return best_map;

  double best_res = std::numeric_limits<double>::infinity();
  int best_shift = -1;

  auto evaluateCandidate = [&](int shift) {
    std::vector<int> test_map = best_map;
    if (!test_map.empty()) test_map[0] = 0;

    const int N = (int)src_fan.size();
    for (int i = 0; i < N; i++)
    {
      const int s_ball = src_fan[i];
      const int j = (shift + i) % N;
      const int d_ball = dst_fan[j];

      if (s_ball < 0 || s_ball >= (int)src.set_to_type.size()) return;
      if (d_ball < 0 || d_ball >= (int)dst.set_to_type.size()) return;

      const int s_ti = src.set_to_type[s_ball];
      const int d_ti = dst.set_to_type[d_ball];
      if (s_ti < 0 || d_ti < 0) return;
      if (s_ti >= (int)test_map.size() || d_ti >= (int)dst.type_to_set.size()) return;
      test_map[s_ti] = d_ti;
    }

    std::vector<std::pair<int,int>> pairs;
    pairs.reserve(1 + src_fan.size());
    pairs.push_back({src.type_to_set[0], dst.type_to_set[test_map[0]]});
    for (int s_ball : src_fan)
    {
      int s_ti = src.set_to_type[s_ball];
      if (s_ti <= 0 || s_ti >= (int)test_map.size()) continue;
      int d_ti = test_map[s_ti];
      if (d_ti < 0 || d_ti >= (int)dst.type_to_set.size()) continue;
      pairs.push_back({src.type_to_set[s_ti], dst.type_to_set[d_ti]});
    }

    Mat5 M;
    double res = 0.0;
    if (fitMobiusFromPairs(*src.parent_set, *dst.parent_set, pairs, M, res) && res < best_res)
    {
      best_res = res;
      best_map = test_map;
      best_shift = shift;
    }
  };

  for (int shift = 0; shift < (int)dst_fan.size(); shift++)
    evaluateCandidate(shift);

  const bool is_A_to_C = (src.parent_set != nullptr && dst.parent_set != nullptr
                          && src.parent_set->name == "icosahedron"
                          && (dst.parent_set->name == "hill-test" || dst.parent_set->name == "hill-testb")
                          && src.type_to_set.size() > 0 && dst.type_to_set.size() > 0
                          && src.type_to_set[0] == 1);
  if (is_A_to_C)
  {
    printFanOrderingDiagnostics(*src.parent_set, src0, src_fan, "A (icosahedron)");
    printFanOrderingDiagnostics(*dst.parent_set, dst0, dst_fan,
                                "C (" + dst.parent_set->name + ")");

    std::cout << "[best-map A->C] shift=" << best_shift
              << " best_res=" << best_res << "\n";
    std::cout << "[best-map A->C] set-id map:";
    for (int sti = 0; sti < (int)best_map.size(); sti++)
    {
      int dti = best_map[sti];
      if (sti < 0 || sti >= (int)src.type_to_set.size()) continue;
      if (dti < 0 || dti >= (int)dst.type_to_set.size()) continue;
      std::cout << " " << src.type_to_set[sti] << "->" << dst.type_to_set[dti];
    }
    std::cout << "\n";
  }

  return best_map;
}

static void rotationAxisAngleFromR(const Eigen::Matrix3d &R,
                                   Eigen::Vector3d &axis,
                                   double &angle)
{
  double c = std::clamp((R.trace() - 1.0) * 0.5, -1.0, 1.0);
  angle = std::acos(c);

  if (angle < 1e-10)
  {
    axis = Eigen::Vector3d(1.0, 0.0, 0.0);
    return;
  }

  Eigen::Vector3d v(R(2,1) - R(1,2),
                    R(0,2) - R(2,0),
                    R(1,0) - R(0,1));
  if (v.norm() > 1e-10)
  {
    axis = v.normalized();
    return;
  }

  // Near-180-degree fallback.
  axis.x() = std::sqrt(std::max(0.0, (R(0,0) + 1.0) * 0.5));
  axis.y() = std::sqrt(std::max(0.0, (R(1,1) + 1.0) * 0.5));
  axis.z() = std::sqrt(std::max(0.0, (R(2,2) + 1.0) * 0.5));
  if (R(0,1) < 0.0) axis.y() = -axis.y();
  if (R(0,2) < 0.0) axis.z() = -axis.z();
  if (axis.norm() < 1e-12) axis = Eigen::Vector3d(1.0, 0.0, 0.0);
  else axis.normalize();
}

static void printMobiusDecomposition(const std::string &label,
                                     const Landscape::Set::Ball::Mobius &mob)
{
  Eigen::Vector3d axis;
  double angle = 0.0;
  rotationAxisAngleFromR(mob.R, axis, angle);
  std::cout << label
            << " similarity=" << (mob.is_similarity ? 1 : 0)
            << " scale=" << mob.s
            << " angle_rad=" << angle
            << " angle_deg=" << (180.0 * angle / pi)
            << " axis=(" << axis.transpose() << ")"
            << " T=(" << mob.T.transpose() << ")"
            << " |T|=" << mob.T.norm()
            << " C=(" << mob.C.transpose() << ")"
            << " |C|=" << mob.C.norm()
            << "\n";
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
                                   Landscape::Set::Ball &dst,
                                   bool quiet = false,
                                   const std::vector<int> *dst_ti_for_src_ti = nullptr)
{
  std::cout << "compute mobius transform" << std::endl;
  src.mobius.M = Mat5::Identity();

  int m = (int)src.type_to_set.size(); // index 0 = ball itself
  if (dst_ti_for_src_ti)
    m = std::min(m, (int)dst_ti_for_src_ti->size());
  const std::string tag = "[mobius " + src.parent_set->name
                        + " ball " + std::to_string(src.type_to_set[0]) + "]";

  if (m < 2)
  {
    if (!quiet) std::cout << tag << " identity (isolated)\n";
    return true;
  }

  Landscape::Set *src_set = src.parent_set;
  Landscape::Set *dst_set = dst.parent_set;

  Mat5 eta = Mat5::Identity();
  eta(4,4) = -1.0;

  std::vector<int> src_tis, dst_tis;
  src_tis.reserve(m);
  dst_tis.reserve(m);
  for (int i = 0; i < m; i++)
  {
    int dti = dst_ti_for_src_ti ? (*dst_ti_for_src_ti)[i] : i;
    if (i < 0 || i >= (int)src.type_to_set.size()) continue;
    if (dti < 0 || dti >= (int)dst.type_to_set.size()) continue;
    src_tis.push_back(i);
    dst_tis.push_back(dti);
  }
  m = (int)src_tis.size();
  if (m < 2)
  {
    if (!quiet) std::cout << tag << " identity (insufficient mapped pairs)\n";
    return true;
  }

  // Build 5×m matrices A (source) and C (dest), columns = unit σ̂.
  Eigen::MatrixXd A(5, m), C(5, m);
  for (int i = 0; i < m; i++)
  {
    const auto &sb = src_set->balls[src.type_to_set[src_tis[i]]];
    const auto &db = dst_set->balls[dst.type_to_set[dst_tis[i]]];
    A.col(i) = conformal_ball(sb);
    C.col(i) = conformal_ball(db);
  }

  // Gram matrix compatibility check.  Print a warning for mismatches but
  // continue — a small mismatch (solver imprecision) gives an approximate M.
  double gram_err = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      gram_err = std::max(gram_err, std::abs(
          mink_dot(A.col(i), A.col(j)) - mink_dot(C.col(i), C.col(j))));
  if (!quiet && gram_err > 1e-3)
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
          int si = src.type_to_set[src_tis[i]], sj = src.type_to_set[src_tis[j]];
          int di = dst.type_to_set[dst_tis[i]], dj = dst.type_to_set[dst_tis[j]];
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

  const double tol_ok   = 1e-5;
  const double tol_fail = 0.1;

  // Project M onto O(4,1) via damped Schulz iteration so that the stored
  // matrix is a valid Möbius transform even when the linear solve left it
  // with residual O(4,1) error.  The projected matrix Mort is used for both
  // transformPoint and the T,C,s,R decomposition so the two are consistent.
  double oo1_after = oo1_err, decomp_err = 0.0;
  if (M.allFinite()) {
    const Mat5 Mort = eta_orthonormalize(M);
    oo1_after = (Mort.transpose() * eta * Mort - eta).norm();

    src.mobius.M = Mort;
    decomposeMobius(Mort, src.mobius);

    // Verify internal consistency: transformPoint (uses Mort) must match
    // transformDecomposed (uses T,C,s,R extracted from Mort).
    if (!quiet) {
      static const Eigen::Vector3d test_pts[] = {
        {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {0.5,0.3,-0.7}, {2,-1,3}
      };
      for (const auto &pt : test_pts)
      {
        auto ref  = src.mobius.transformPoint(pt);
        auto fast = src.mobius.transformDecomposed(pt);
        if (ref.allFinite() && fast.allFinite())
          decomp_err = std::max(decomp_err, (ref - fast).norm());
      }
    }
  }

  // Treat NaN as FAIL: use oo1_after (post-Schulz) since that is the stored
  // value — the raw oo1_err can be large even when projection succeeds.
  bool bad = !std::isfinite(oo1_after) || !std::isfinite(residual)
             || oo1_after > tol_fail || residual > tol_fail;

  if (!quiet && (residual > tol_ok || bad))
  {
    std::cout << tag << "  r=" << r
              << "  gram_err=" << gram_err
              << "  O(4,1)_err=" << oo1_err
              << "  oo1_after=" << oo1_after
              << "  residual=" << residual
              << "  decomp_err=" << decomp_err;
  }

  if (bad)
  {
    if (!quiet) std::cout << "  FAIL\n";
    return false;
  }

  if (!quiet && residual > tol_ok)
    std::cout << ((oo1_after > tol_ok || residual > tol_ok) ? "  APPROX\n" : "  OK\n");
  return true;
}


void Landscape::printConnectivity(bool show_valid_destinations)
{
  for (auto &set : sets)
  {
    std::cout << set.name << ": " << set.balls.size() << " balls" << std::endl;
    for (int i = 0; i<(int)set.balls.size(); i++)
    {
      auto &ball = set.balls[i];

      if (show_valid_destinations)
      {
        bool written_ball = false;
        Set *last_par = nullptr;
        for (auto *dest_ball: ball.type->balls)
        {
          Set *par = dest_ball->parent_set;
          if (par != ball.parent_set)
          {
            for (int j = 0; j<par->balls.size(); j++)
            {
              if (&par->balls[j] == dest_ball)
              {
                if (!written_ball)
                {
                  std::cout << " ball " << i << ":";
                  written_ball = true;
                }
                if (par != last_par)
                  std::cout << " " << par->name << ": ";
                std::cout << j << ",";
              }
            }
          }
          last_par = par;
        }
        if (written_ball)
          std::cout << std::endl;
      }
      Set::Ball *dest_ball = ball.dest_ball;
      if (dest_ball != nullptr && dest_ball != &ball)
      {
        if (!show_valid_destinations)
          std::cout << " ball " << i << ":";
        Set *par = dest_ball->parent_set;
        std::cout << "  connects to " << par->name;
        if (ball.dest_ball_id == -1)
          std::cout << " first valid";
        std::cout << " ball: ";
        for (int j = 0; j<par->balls.size(); j++)
          if (&par->balls[j] == dest_ball)
            std::cout << j;
        std::cout << std::endl;
      }
    }
  }
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
  // Step 1: resolve dest_ball pointers by name.
  for (auto &set : sets)
  {
    for (auto &ball : set.balls)
    {
      if (ball.dest_set == "")
      {
        ball.dest_ball = &ball;
        continue;
      }

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
        ball.dest_ball = &dest_set_ptr->balls[0];
    }
  }

  const double damping = 1e-10;
  const double k = 0.0;
  const double sor = 1.2;
  const double anchor_weight = 0.0;
  const bool fit_mobius_transform = true;

  std::mt19937 rng(42);
  bool warmup_done = !fit_mobius_transform;

  // Within-set pairs.
  struct IntraPair { int si, i, j; };
  std::vector<IntraPair> intra_pairs;
  std::set<std::tuple<int,int,int>> intra_pair_keys;
  auto addIntraPair = [&](int si, int a, int b) {
    if (si < 0 || si >= (int)sets.size()) return;
    if (a < 0 || b < 0 || a >= (int)sets[si].balls.size() || b >= (int)sets[si].balls.size()) return;
    int i = std::max(a, b), j = std::min(a, b);
    auto key = std::make_tuple(si, i, j);
    if (intra_pair_keys.insert(key).second)
      intra_pairs.push_back({si, i, j});
  };
  for (int si = 0; si < (int)sets.size(); si++)
  {
    int n = (int)sets[si].balls.size();
    for (int i = 0; i < n; i++)
      for (int j = 0; j < i; j++)
        addIntraPair(si, i, j);
  }

  auto setIndexByName = [&](const std::string &name) -> int {
    for (int si = 0; si < (int)sets.size(); si++)
      if (sets[si].name == name) return si;
    return -1;
  };

  // Cross-set Möbius constraints.
  struct MobiusPair { Set::Ball *src; Set::Ball *dst; int src_ti; int dst_ti; };
  struct MobiusLink { Set::Ball *src; Set::Ball *dst; std::vector<int> dst_ti_for_src_ti; };
  struct OverlapLink {
    Set::Ball *ov_ball;
    Set::Ball *dst_ball;
    Set::Ball *mobius_owner;
    int src_anchor_set_ball;
    int dst_anchor_set_ball;
  };

  std::vector<MobiusPair> mobius_pairs;
  std::vector<MobiusLink> mobius_links;
  std::set<std::tuple<Set::Ball*, Set::Ball*>> mobius_link_keys;
  std::vector<OverlapLink> overlap_links;

  auto addMobiusLink = [&](Set::Ball *src, Set::Ball *dst) {
    if (src == nullptr || dst == nullptr || src == dst) return;
    auto key = std::make_tuple(src, dst);
    if (mobius_link_keys.insert(key).second)
      mobius_links.push_back({src, dst, identityTypeMap(*src, *dst)});
  };

  for (auto &set : sets)
    for (auto &ball : set.balls)
      if (ball.dest_ball != nullptr && ball.dest_ball != &ball)
        addMobiusLink(&ball, ball.dest_ball);

  auto findLinkMap = [&](Set::Ball *src, Set::Ball *dst) -> const std::vector<int>* {
    for (auto &lk : mobius_links)
      if (lk.src == src && lk.dst == dst)
        return &lk.dst_ti_for_src_ti;
    return nullptr;
  };

  auto buildOverlapLinks = [&]() {
    overlap_links.clear();
    for (int si = 0; si < (int)sets.size(); si++)
    {
      Set &src_set = sets[si];
      for (const auto &ov : src_set.overlaps)
      {
        if (ov.ball_0 < 0 || ov.ball_0 >= (int)src_set.balls.size()
            || ov.ball_1 < 0 || ov.ball_1 >= (int)src_set.balls.size())
          continue;

        Set::Ball &parent0 = src_set.balls[ov.ball_0];
        Set::Ball &parent1 = src_set.balls[ov.ball_1];

        if (parent0.dest_ball == nullptr || parent1.dest_ball == nullptr) continue;
        if (parent0.dest_ball == &parent0 || parent1.dest_ball == &parent1) continue;

        const std::vector<int> *map0 = findLinkMap(&parent0, parent0.dest_ball);
        const std::vector<int> *map1 = findLinkMap(&parent1, parent1.dest_ball);
        if (map0 == nullptr || map1 == nullptr) continue;

        // For D[a], constrain against the equivalent of A[a] on parent1's dest set,
        // using the A[b]→dest match as orientation anchor.
        int dst_for_a = mapSetBallThroughTypeMap(parent1, *parent1.dest_ball, *map1, ov.ball_0);
        int dst_for_b = mapSetBallThroughTypeMap(parent0, *parent0.dest_ball, *map0, ov.ball_1);
        if (dst_for_a < 0 || dst_for_b < 0) continue;

        Set *dst_set_for_a = parent1.dest_ball->parent_set;
        Set *dst_set_for_b = parent0.dest_ball->parent_set;
        if (dst_set_for_a == nullptr || dst_set_for_b == nullptr) continue;

        if (dst_for_a >= (int)dst_set_for_a->balls.size()) continue;
        if (dst_for_b >= (int)dst_set_for_b->balls.size()) continue;

        int dsiO = setIndexByName(ov.dest_set);
        if (dsiO < 0) continue;

        Set &ov_set = sets[dsiO];
        Set &dstA = *dst_set_for_a;
        Set &dstB = *dst_set_for_b;

        int a = ov.ball_0;
        int b = ov.ball_1;
        if (a < 0 || a >= (int)ov_set.balls.size() || b < 0 || b >= (int)ov_set.balls.size()) continue;

        const int anchor_a_dst = parent1.dest_ball->type_to_set[0];
        const int anchor_b_dst = parent0.dest_ball->type_to_set[0];

        // D[a] is constrained to parent1-dest equivalent of A[a], anchored by A[b].
        overlap_links.push_back({&ov_set.balls[a], &dstA.balls[dst_for_a], &parent1,
                                 /*src_anchor_set_ball=*/b,
                                 /*dst_anchor_set_ball=*/anchor_a_dst});
        // D[b] is constrained to parent0-dest equivalent of A[b], anchored by A[a].
        overlap_links.push_back({&ov_set.balls[b], &dstB.balls[dst_for_b], &parent0,
                                 /*src_anchor_set_ball=*/a,
                                 /*dst_anchor_set_ball=*/anchor_b_dst});
      }
    }
  };

  auto rebuildMobiusPairs = [&]() {
    mobius_pairs.clear();
    std::set<std::tuple<Set::Ball*, Set::Ball*, int, int>> mobius_pair_keys;
    for (auto &lk : mobius_links)
    {
      int m = (int)lk.dst_ti_for_src_ti.size();
      for (int sti = 0; sti < m; sti++)
      {
        int dti = lk.dst_ti_for_src_ti[sti];
        if (sti < 0 || sti >= (int)lk.src->type_to_set.size()) continue;
        if (dti < 0 || dti >= (int)lk.dst->type_to_set.size()) continue;
        auto key = std::make_tuple(lk.src, lk.dst, sti, dti);
        if (mobius_pair_keys.insert(key).second)
          mobius_pairs.push_back({lk.src, lk.dst, sti, dti});
      }
    }
  };
  rebuildMobiusPairs();

  if (!fit_mobius_transform)
  {
    for (auto &lk : mobius_links)
      lk.src->mobius = Set::Ball::Mobius();
  }

  std::vector<std::vector<Eigen::Vector3d>> ref_centres(sets.size());
  std::vector<std::vector<double>> ref_radii(sets.size());
  for (int si = 0; si < (int)sets.size(); si++)
  {
    int n = (int)sets[si].balls.size();
    ref_centres[si].resize(n);
    ref_radii[si].resize(n);
    for (int i = 0; i < n; i++)
    {
      ref_centres[si][i] = sets[si].balls[i].centre;
      ref_radii[si][i] = clampSignedRadius(sets[si].balls[i].radius);
    }
  }

  auto sigma_to_sphere_grad = [](double r, const Eigen::Vector3d &C, const Vec5 &gs)
      -> std::pair<Eigen::Vector3d, double>
  {
    double C2 = C.squaredNorm();
    double rs = clampSignedRadius(r);
    Eigen::Vector3d gC = (gs.head<3>() + (gs(4) - gs(3)) * C) / rs;
    double gr = (-C / (rs*rs)).dot(gs.head<3>())
              + (rs*rs - 1.0 + C2) / (2.0*rs*rs) * gs(3)
              - (rs*rs + 1.0 + C2) / (2.0*rs*rs) * gs(4);
    return {gC, gr};
  };

  auto reportWarmupDiagnostics = [&]() {
    std::cout << "[warmup-end] ===== constraint diagnostics =====\n";

    struct Stat {
      int n = 0;
      double sum2 = 0.0;
      double max_abs = 0.0;
      void add(double v)
      {
        double a = std::abs(v);
        n++;
        sum2 += v * v;
        if (a > max_abs) max_abs = a;
      }
    };

    auto rms = [](const Stat &s) {
      return s.n > 0 ? std::sqrt(s.sum2 / (double)s.n) : 0.0;
    };

    // 1) Intra-constraint errors by set.
    std::vector<Stat> intra_stats(sets.size());
    for (int si = 0; si < (int)sets.size(); si++)
    {
      auto &set = sets[si];
      int n = (int)set.balls.size();
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < i; j++)
        {
          int order = set.conn(i, j);
          const auto &bi = set.balls[i];
          const auto &bj = set.balls[j];

          double ri = clampSignedRadius(bi.radius);
          double rj = clampSignedRadius(bj.radius);
          double d = (bi.centre - bj.centre).norm();
          double residual = 0.0;

          if (order <= 0)
          {
            // Same quantity as solver branch: err = d - (ri+rj+k_{order==0}).
            double err = d - (ri + rj + (order == 0 ? k : 0.0));
            residual = (order == 0) ? std::min(0.0, err) : err;
          }
          else
          {
            double cos_theta = (d*d - ri*ri - rj*rj) / (2.0 * ri * rj);
            double theta = std::acos(std::clamp(cos_theta, -1.0, 1.0));
            residual = pi / (double)order - theta;
          }

          intra_stats[si].add(residual);
        }
      }
    }

    for (int si = 0; si < (int)sets.size(); si++)
    {
      std::cout << "[warmup-end][intra-summary] " << sets[si].name
                << " n=" << intra_stats[si].n
                << " rms=" << rms(intra_stats[si])
                << " max_abs=" << intra_stats[si].max_abs << "\n";
    }

    // Classify overlap links so we can separate primary inter constraints from cross-overlap constraints.
    std::set<std::tuple<Set::Ball*, Set::Ball*>> overlap_link_keys;
    for (const auto &ol : overlap_links)
      overlap_link_keys.insert(std::make_tuple(ol.ov_ball, ol.dst_ball));

    Stat inter_stats, cross_stats;
    std::map<std::string, Stat> inter_link_stats;
    std::map<std::string, Stat> cross_link_stats;

    // 2) Inter-constraints and 3) overlap cross-constraints.
    for (const auto &mp : mobius_pairs)
    {
      Set::Ball &src_ball = *mp.src;
      Set::Ball &dst_ball = *mp.dst;

      int si = src_ball.type_to_set[mp.src_ti];
      int di = dst_ball.type_to_set[mp.dst_ti];
      Set::Ball &sA = src_ball.parent_set->balls[si];
      Set::Ball &dA = dst_ball.parent_set->balls[di];

      bool is_cross = overlap_link_keys.count(std::make_tuple(mp.src, mp.dst)) > 0;

      const Set::Ball::Mobius *mobius_to_use = &src_ball.mobius;
      if (is_cross)
      {
        for (const auto &ol : overlap_links)
        {
          if (ol.ov_ball == mp.src && ol.dst_ball == mp.dst)
          {
            mobius_to_use = &ol.mobius_owner->mobius;
            break;
          }
        }
      }

      Vec5 sigma_s = conformal_sphere(sA.centre, clampSignedRadius(sA.radius));
      Vec5 sigma_d = conformal_sphere(dA.centre, clampSignedRadius(dA.radius));
      double residual = (mobius_to_use->M * sigma_s - sigma_d).norm();

      std::string link_key = src_ball.parent_set->name + " ball " + std::to_string(src_ball.type_to_set[0])
                           + " -> " + dst_ball.parent_set->name + " ball " + std::to_string(dst_ball.type_to_set[0]);

      if (is_cross)
      {
        cross_stats.add(residual);
        cross_link_stats[link_key].add(residual);
      }
      else
      {
        inter_stats.add(residual);
        inter_link_stats[link_key].add(residual);
      }
    }

    for (const auto &kv : inter_link_stats)
    {
      std::cout << "[warmup-end][inter-summary] " << kv.first
                << " n=" << kv.second.n
                << " rms=" << rms(kv.second)
                << " max_abs=" << kv.second.max_abs << "\n";
    }

    for (const auto &kv : cross_link_stats)
    {
      std::cout << "[warmup-end][cross-summary] " << kv.first
                << " n=" << kv.second.n
                << " rms=" << rms(kv.second)
                << " max_abs=" << kv.second.max_abs << "\n";
    }

    std::cout << "[warmup-end][inter-total] n=" << inter_stats.n
              << " rms=" << rms(inter_stats)
              << " max_abs=" << inter_stats.max_abs << "\n";
    std::cout << "[warmup-end][cross-total] n=" << cross_stats.n
              << " rms=" << rms(cross_stats)
              << " max_abs=" << cross_stats.max_abs << "\n";
    std::cout << "[warmup-end] =================================\n";
  };

  const int warmup_iters = 1500;
  const int effective_warmup_iters = fit_mobius_transform ? warmup_iters : 0;
  for (int it = 0; it < iterations + effective_warmup_iters; it++)
  {
    if (fit_mobius_transform && !warmup_done && it >= warmup_iters)
    {
      for (auto &lk : mobius_links)
      {
        lk.dst_ti_for_src_ti = findBestNeighbourTypeMap(*lk.src, *lk.dst);
        computeMobiusTransform(*lk.src, *lk.dst, /*quiet=*/true, &lk.dst_ti_for_src_ti);

        const bool is_A_to_C = (lk.src != nullptr && lk.dst != nullptr
                                && lk.src->parent_set != nullptr && lk.dst->parent_set != nullptr
                                && lk.src->parent_set->name == "icosahedron"
              && (lk.dst->parent_set->name == "hill-test" || lk.dst->parent_set->name == "hill-testb")
                                && lk.src->type_to_set.size() > 0 && lk.dst->type_to_set.size() > 0
              && lk.src->type_to_set[0] == 1);
        if (is_A_to_C)
        {
          printMobiusDecomposition("[mobius A->C]", lk.src->mobius);
        }
      }

      // Add overlap links only after warmup so fitting is unaffected by overlaps.
      buildOverlapLinks();

      for (const auto &ol : overlap_links)
      {
        addMobiusLink(ol.ov_ball, ol.dst_ball);

        // Overlap links must inherit neighbour mapping from the corresponding
        // parent cross-set link (ball_0's link-set or ball_1's link-set),
        // rather than running an independent best-fit search.
        const std::vector<int> *owner_map = nullptr;
        if (ol.mobius_owner != nullptr
            && ol.mobius_owner->dest_ball != nullptr
            && ol.mobius_owner->dest_ball != ol.mobius_owner)
        {
          owner_map = findLinkMap(ol.mobius_owner, ol.mobius_owner->dest_ball);
        }

        for (auto &lk : mobius_links)
        {
          if (lk.src == ol.ov_ball && lk.dst == ol.dst_ball)
          {
            if (owner_map != nullptr)
            {
              // Overlap spheres D are in A-coordinates, but each D->dst link
              // must follow the chosen cyclic shift for that branch (A->B/A->C).
              // Build the local fan map from the selected dst center first,
              // anchored by the paired overlap parent, then only fill gaps from
              // owner inheritance.
              lk.dst_ti_for_src_ti = buildFanTypeMap(*lk.src, *lk.dst,
                                                     ol.src_anchor_set_ball,
                                                     ol.dst_anchor_set_ball);

              std::vector<int> inherited = buildInheritedTypeMapBySetIndex(
                  *lk.src, *lk.dst, *ol.mobius_owner, *owner_map);
              for (int ti = 0; ti < (int)lk.dst_ti_for_src_ti.size() && ti < (int)inherited.size(); ti++)
              {
                if (lk.dst_ti_for_src_ti[ti] < 0 && inherited[ti] >= 0)
                  lk.dst_ti_for_src_ti[ti] = inherited[ti];
              }
            }
            else
              lk.dst_ti_for_src_ti = identityTypeMap(*lk.src, *lk.dst);

            if (lk.src != nullptr && lk.dst != nullptr
                && lk.src->parent_set != nullptr && lk.dst->parent_set != nullptr
                && lk.src->parent_set->name == "tree-hill-test"
                && !lk.src->type_to_set.empty())
            {
              const int d_center = lk.src->type_to_set[0];
              if (d_center == 0 || d_center == 1)
              {
                std::string label = "D(" + std::to_string(d_center) + ")->"
                                  + lk.dst->parent_set->name
                                  + "(" + std::to_string(lk.dst->type_to_set[0]) + ")";
                printMappedFanAdjacencyDiagnostics(*lk.src, *lk.dst,
                                                   lk.dst_ti_for_src_ti,
                                                   label);
              }
            }
            break;
          }
        }
      }

      // Targeted debug: report overlap-link destination ids for requested balls.
      auto printMappedOverlapBall = [&](const std::string &src_set_name,
                                        int src_ball_id,
                                        const std::string &dst_set_name)
      {
        for (const auto &lk : mobius_links)
        {
          if (lk.src == nullptr || lk.dst == nullptr) continue;
          if (lk.src->parent_set == nullptr || lk.dst->parent_set == nullptr) continue;
          if (lk.src->parent_set->name != src_set_name) continue;
          if (lk.dst->parent_set->name != dst_set_name) continue;

          int mapped = mapSetBallThroughTypeMap(*lk.src, *lk.dst, lk.dst_ti_for_src_ti, src_ball_id);
          if (mapped >= 0)
          {
            std::cout << "[overlap-map] " << src_set_name << " ball " << src_ball_id
                      << " links to " << dst_set_name << " ball " << mapped
                      << " (via src center ball " << lk.src->type_to_set[0]
                      << " -> dst center ball " << lk.dst->type_to_set[0] << ")\n";
            return;
          }
        }
        std::cout << "[overlap-map] no mapping found for " << src_set_name
                  << " ball " << src_ball_id << " to " << dst_set_name << "\n";
      };

      printMappedOverlapBall("tree-hill-test", 6, "tree-test");
      printMappedOverlapBall("tree-hill-test", 3, "hill-testb");

      rebuildMobiusPairs();
      reportWarmupDiagnostics();
      warmup_done = true;
    }

    // Keep warmup deterministic and per-topology stable: this prevents
    // unrelated additions (e.g. overlap-only sets) from changing the
    // A->C fit simply by perturbing global shuffle order.
    if (warmup_done)
      std::shuffle(intra_pairs.begin(), intra_pairs.end(), rng);
    std::shuffle(mobius_pairs.begin(), mobius_pairs.end(), rng);

    // (A) intra-set constraints
    for (auto [si, i, j] : intra_pairs)
    {
      Set &set = sets[si];
      int order = set.conn(i, j);
      Set::Ball &bi = set.balls[i];
      Set::Ball &bj = set.balls[j];

      Eigen::Vector3d Ci = bi.centre;
      Eigen::Vector3d Cj = bj.centre;
      double ri = clampSignedRadius(bi.radius);
      double rj = clampSignedRadius(bj.radius);

      Eigen::Vector3d Delta = Ci - Cj;
      double d = Delta.norm();
      if (d < 1e-12) continue;

      double err = 0.0;
      Eigen::Vector3d gCi = Eigen::Vector3d::Zero(), gCj = Eigen::Vector3d::Zero();
      double gri = 0.0, grj = 0.0;

      if (order <= 0)
      {
        double targ_d = ri + rj + (order == 0 ? k : 0.0);
        err = d - targ_d;
        if (order == 0 && err >= 0.0) continue;
        gCi =  Delta / d;
        gCj = -Delta / d;
        gri = -1.0;
        grj = -1.0;
      }
      else
      {
        double d2 = d * d;
        double cos_theta = (d2 - ri*ri - rj*rj) / (2.0 * ri * rj);
        if (cos_theta >= 1.0 || cos_theta <= -1.0)
        {
          err = d - (ri + rj);
          gCi =  Delta / d;
          gCj = -Delta / d;
          gri = -1.0;
          grj = -1.0;
        }
        else
        {
          double theta = std::acos(cos_theta);
          double sin_theta = std::sin(theta);
          if (std::abs(sin_theta) < 1e-10) continue;
          err = pi / (double)order - theta;
          double inv_sin = 1.0 / sin_theta;
          gCi =  inv_sin * Delta / (ri * rj);
          gCj = -inv_sin * Delta / (ri * rj);
          gri = inv_sin * (-(ri*ri + d2 - rj*rj) / (2.0 * ri*ri * rj));
          grj = inv_sin * (-(rj*rj + d2 - ri*ri) / (2.0 * rj*rj * ri));
        }
      }

      const double wi = bi.mobility;
      const double wj = bj.mobility;
      gCi *= wi;  gri *= wi;
      gCj *= wj;  grj *= wj;

      double g2 = gCi.squaredNorm() + gri*gri + gCj.squaredNorm() + grj*grj;
      double step = -sor * err / (g2 + damping);

      bi.centre = Ci + step * gCi;
      bi.radius = clampSignedRadius(ri + step * gri);
      bj.centre = Cj + step * gCj;
      bj.radius = clampSignedRadius(rj + step * grj);
    }

    // (C) proximity anchors (apply during warmup too, to fix gauge drift)
    const double anchor_step = sor * anchor_weight;
    for (int si = 0; si < (int)sets.size(); si++)
    {
      int n = (int)sets[si].balls.size();
      for (int i = 0; i < n; i++)
      {
        Set::Ball &b = sets[si].balls[i];
        const Eigen::Vector3d &C0 = ref_centres[si][i];
        double r0 = ref_radii[si][i];

        double w = b.mobility;
        if (w <= 0.0) continue;

        b.centre -= (anchor_step * w) * (b.centre - C0);
        b.radius = clampSignedRadius(b.radius - (anchor_step * w) * (b.radius - r0));
      }
    }

    if (!warmup_done) continue;

    // (B) cross-set Möbius residual constraints
    for (auto &mp : mobius_pairs)
    {
      Set::Ball &src_ball = *mp.src;
      Set::Ball &dst_ball = *mp.dst;

      int si = src_ball.type_to_set[mp.src_ti];
      int di = dst_ball.type_to_set[mp.dst_ti];

      Set::Ball &sA = src_ball.parent_set->balls[si];
      Set::Ball &dA = dst_ball.parent_set->balls[di];

      const Set::Ball::Mobius *mobius_to_use = &src_ball.mobius;
      for (const auto &ol : overlap_links)
      {
        if (ol.ov_ball == mp.src && ol.dst_ball == mp.dst)
        {
          mobius_to_use = &ol.mobius_owner->mobius;
          break;
        }
      }

      double rS = clampSignedRadius(sA.radius), rD = clampSignedRadius(dA.radius);
      Eigen::Vector3d CS = sA.centre, CD = dA.centre;

      Vec5 sigma_s = conformal_sphere(CS, rS);
      Vec5 sigma_d = conformal_sphere(CD, rD);

      Vec5 e = mobius_to_use->M * sigma_s - sigma_d;
      Vec5 g_sigma_s = mobius_to_use->M.transpose() * e;
      Vec5 g_sigma_d = -e;

      auto [gCS, gRS] = sigma_to_sphere_grad(rS, CS, g_sigma_s);
      auto [gCD, gRD] = sigma_to_sphere_grad(rD, CD, g_sigma_d);

      const double wS = sA.mobility;
      const double wD = dA.mobility;
      gCS *= wS;  gRS *= wS;
      gCD *= wD;  gRD *= wD;

      double g2 = gCS.squaredNorm() + gRS*gRS + gCD.squaredNorm() + gRD*gRD;
      double step = -sor * e.squaredNorm() / (2.0 * (g2 + damping));

      sA.centre += step * gCS;
      sA.radius = clampSignedRadius(sA.radius + step * gRS);
      dA.centre += step * gCD;
      dA.radius = clampSignedRadius(dA.radius + step * gRD);
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
        // computeMobiusTransform(ball, ball); // no need to compute as we won't be using the transform if no dest
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
