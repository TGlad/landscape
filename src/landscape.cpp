#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

static const double pi = std::acos(-1.0);

#include <iostream>
#include <vector>
#include <set>
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

static std::vector<int> findBestNeighbourTypeMap(const Landscape::Set::Ball &src,
                                                 const Landscape::Set::Ball &dst)
{
  std::vector<int> best_map = identityTypeMap(src, dst);
  int m = (int)best_map.size();
  if (m < 2) return best_map;

  int src0 = src.type_to_set[0];
  int dst0 = dst.type_to_set[0];

  std::vector<int> src_neigh_tis, dst_neigh_tis;
  for (int ti = 1; ti < m; ti++)
  {
    if (src.parent_set->conn(src0, src.type_to_set[ti]) > 0)
      src_neigh_tis.push_back(ti);
    if (dst.parent_set->conn(dst0, dst.type_to_set[ti]) > 0)
      dst_neigh_tis.push_back(ti);
  }

  if (src_neigh_tis.empty() || src_neigh_tis.size() != dst_neigh_tis.size())
    return best_map;

  auto sideFilterPasses = [&](const std::vector<int> &test_map) -> bool {
    if (src_neigh_tis.size() < 3) return true;

    const Landscape::Set &src_set = *src.parent_set;
    const Landscape::Set &dst_set = *dst.parent_set;

    std::vector<Eigen::Vector3d> src_pts, dst_pts;
    src_pts.reserve(src_neigh_tis.size());
    dst_pts.reserve(src_neigh_tis.size());
    for (int s_ti : src_neigh_tis)
    {
      int d_ti = test_map[s_ti];
      src_pts.push_back(src_set.balls[src.type_to_set[s_ti]].centre);
      dst_pts.push_back(dst_set.balls[dst.type_to_set[d_ti]].centre);
    }

    auto meanPoint = [](const std::vector<Eigen::Vector3d> &pts) -> Eigen::Vector3d {
      Eigen::Vector3d c = Eigen::Vector3d::Zero();
      for (const auto &p : pts) c += p;
      return c / (double)std::max(1, (int)pts.size());
    };

    auto setCentroid = [](const Landscape::Set &s) -> Eigen::Vector3d {
      if (s.balls.empty()) return Eigen::Vector3d::Zero();
      Eigen::Vector3d c = Eigen::Vector3d::Zero();
      for (const auto &b : s.balls) c += b.centre;
      return c / (double)s.balls.size();
    };

    auto orientedPlaneNormal = [&](const std::vector<Eigen::Vector3d> &pts,
                                   const Eigen::Vector3d &center_ball,
                                   Eigen::Vector3d &origin,
                                   Eigen::Vector3d &normal) -> bool {
      origin = meanPoint(pts);
      Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
      for (const auto &p : pts)
      {
        Eigen::Vector3d d = p - origin;
        cov += d * d.transpose();
      }
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
      if (eig.info() != Eigen::Success) return false;
      normal = eig.eigenvectors().col(0).normalized(); // smallest variance axis
      if (!normal.allFinite() || normal.squaredNorm() < 1e-20) return false;

      // Fix normal sign consistently using center ball side.
      if ((center_ball - origin).dot(normal) < 0.0)
        normal = -normal;
      return true;
    };

    Eigen::Vector3d src_o, src_n, dst_o, dst_n;
    const Eigen::Vector3d src_center_ball = src_set.balls[src.type_to_set[0]].centre;
    const Eigen::Vector3d dst_center_ball = dst_set.balls[dst.type_to_set[0]].centre;
    if (!orientedPlaneNormal(src_pts, src_center_ball, src_o, src_n)) return true;
    if (!orientedPlaneNormal(dst_pts, dst_center_ball, dst_o, dst_n)) return true;

    double src_side = (setCentroid(src_set) - src_o).dot(src_n);
    double dst_side = (setCentroid(dst_set) - dst_o).dot(dst_n);

    // Near-coplanar centroid wrt fan plane: treat as ambiguous and allow.
    const double eps = 1e-8;
    if (std::abs(src_side) < eps || std::abs(dst_side) < eps)
      return true;

    return src_side * dst_side > 0.0;
  };

  std::vector<int> chosen(src_neigh_tis.size(), -1);
  std::vector<char> used(dst_neigh_tis.size(), 0);
  double best_res = std::numeric_limits<double>::infinity();

  std::function<void(int)> dfs = [&](int idx)
  {
    if (idx == (int)src_neigh_tis.size())
    {
      std::vector<int> test_map = best_map;
      for (int k = 0; k < (int)src_neigh_tis.size(); k++)
        test_map[src_neigh_tis[k]] = chosen[k];

      if (!sideFilterPasses(test_map))
        return;

      std::vector<std::pair<int,int>> pairs;
      pairs.reserve(1 + src_neigh_tis.size());
      pairs.push_back({src.type_to_set[0], dst.type_to_set[test_map[0]]});
      for (int s_ti : src_neigh_tis)
        pairs.push_back({src.type_to_set[s_ti], dst.type_to_set[test_map[s_ti]]});

      Mat5 M;
      double res = 0.0;
      if (fitMobiusFromPairs(*src.parent_set, *dst.parent_set, pairs, M, res) && res < best_res)
      {
        best_res = res;
        best_map = test_map;
      }
      return;
    }

    int s_ti = src_neigh_tis[idx];
    for (int c = 0; c < (int)dst_neigh_tis.size(); c++)
    {
      if (used[c]) continue;
      int d_ti = dst_neigh_tis[c];

      bool ok = true;
      for (int p = 0; p < idx; p++)
      {
        int s_prev = src_neigh_tis[p];
        int d_prev = chosen[p];
        int so = src.parent_set->conn(src.type_to_set[s_ti], src.type_to_set[s_prev]);
        int do_ = dst.parent_set->conn(dst.type_to_set[d_ti], dst.type_to_set[d_prev]);
        if (so != do_)
        {
          ok = false;
          break;
        }
      }
      if (!ok) continue;

      chosen[idx] = d_ti;
      used[c] = 1;
      dfs(idx + 1);
      used[c] = 0;
      chosen[idx] = -1;
    }
  };

  dfs(0);
  return best_map;
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

  // Resolve overlap destination sets and collect additional overlap-driven
  // src→dst Möbius links. For overlap (a,b), with
  //   A = balls[a].dest_set, B = balls[b].dest_set,
  // use cross-links:
  //   a -> B[a],  b -> A[b]
  // and neighbour links (same-index assumption):
  //   n neighbour of a (n != b) -> B[n]
  //   n neighbour of b (n != a) -> A[n]
  // If n is neighbour of both a and b, both links are added (dedup later).
  std::vector<std::pair<Set::Ball*, Set::Ball*>> overlap_links;

  // ── Step 2: build all constraint pairs ──────────────────────────────────
  //
  // Each entry is one scalar constraint.  We use a tagged union approach:
  // either a within-set pair (conn order), or a cross-set Gram-matching pair.
  const double damping = 1e-10;
  const double k = 0.0; // minimum extra gap for separation constraints
  const double sor = 1.2; // successive over-relaxation factor
  const double anchor_weight = 0.0;//1e-3; // formal proximity anchor weight

  // Within-set pairs: {set index, ball i, ball j}
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

  // Helpers used by overlap constraint wiring.
  auto setIndexByName = [&](const std::string &name) -> int {
    for (int si = 0; si < (int)sets.size(); si++)
      if (sets[si].name == name) return si;
    return -1;
  };

  for (int si = 0; si < (int)sets.size(); si++)
  {
    Set &src_set = sets[si];
    for (const auto &ov : src_set.overlaps)
    {
      if (ov.ball_0 < 0 || ov.ball_0 >= (int)src_set.balls.size()
          || ov.ball_1 < 0 || ov.ball_1 >= (int)src_set.balls.size())
      {
        std::cerr << "[applyConnectivity] overlap has invalid source ball indices in set '"
                  << src_set.name << "'\n";
        continue;
      }

      Set::Ball &parent0 = src_set.balls[ov.ball_0];
      Set::Ball &parent1 = src_set.balls[ov.ball_1];

      int dsiO = setIndexByName(ov.dest_set);      // overlap set (source of overlap links)
      // If a parent ball has no destination set, fall back to the parent/source set.
      int dsiA = parent0.dest_set.empty() ? si : setIndexByName(parent0.dest_set);
      int dsiB = parent1.dest_set.empty() ? si : setIndexByName(parent1.dest_set);
      if (dsiO < 0)
      {
        std::cerr << "[applyConnectivity] overlap source set '" << ov.dest_set
                  << "' not found for parent set '" << src_set.name << "'\n";
        continue;
      }
      if (dsiA < 0 || dsiB < 0)
      {
        std::cerr << "[applyConnectivity] overlap needs valid destination (or source-fallback) sets on both source balls in set '"
                  << src_set.name << "'\n";
        continue;
      }
      Set &ov_set = sets[dsiO];
      Set &dstA = sets[dsiA];
      Set &dstB = sets[dsiB];

      int a = ov.ball_0;
      int b = ov.ball_1;

      if (a < 0 || a >= (int)ov_set.balls.size() || b < 0 || b >= (int)ov_set.balls.size())
      {
        std::cerr << "[applyConnectivity] overlap indices out of range for overlap set '"
                  << ov_set.name << "'\n";
        continue;
      }

      if (a < 0 || a >= (int)dstB.balls.size())
      {
        std::cerr << "[applyConnectivity] overlap cross-link a->B out of range for set '"
                  << dstB.name << "'\n";
        continue;
      }
      if (b < 0 || b >= (int)dstA.balls.size())
      {
        std::cerr << "[applyConnectivity] overlap cross-link b->A out of range for set '"
                  << dstA.name << "'\n";
        continue;
      }

      // Cross links: a -> B[a], b -> A[b]. Annotate with which parent mobius to use (0 or 1).
      // We use a struct to carry the mobius parent index.
      struct OverlapLink { Set::Ball *ov_ball; Set::Ball *dst_ball; int parent_mobius; };
      static std::vector<OverlapLink> overlap_links_ext; // static to avoid redefinition warning
      overlap_links_ext.clear();
      overlap_links_ext.push_back({&ov_set.balls[a], &dstB.balls[a], 0}); // use parent0.mobius
      overlap_links_ext.push_back({&ov_set.balls[b], &dstA.balls[b], 1}); // use parent1.mobius

      // Rule 3: neighbours of a map to corresponding neighbours of B[a], use parent0.mobius
      for (int n = 0; n < (int)ov_set.balls.size(); n++)
      {
        if (n == b) continue; // handled by anchor link a -> B[a]
        if (ov_set.conn(a, n) > 0)
        {
          if (n >= 0 && n < (int)dstB.balls.size())
            overlap_links_ext.push_back({&ov_set.balls[n], &dstB.balls[n], 0});
          else
            std::cerr << "[applyConnectivity] overlap neighbour index " << n
                      << " out of range for destination set '" << dstB.name << "'\n";
        }
      }

      // Rule 4: neighbours of b map to corresponding neighbours of A[b], use parent1.mobius
      for (int n = 0; n < (int)ov_set.balls.size(); n++)
      {
        if (n == a) continue; // handled by anchor link b -> A[b]
        if (ov_set.conn(b, n) > 0)
        {
          if (n >= 0 && n < (int)dstA.balls.size())
            overlap_links_ext.push_back({&ov_set.balls[n], &dstA.balls[n], 1});
          else
            std::cerr << "[applyConnectivity] overlap neighbour index " << n
                      << " out of range for destination set '" << dstA.name << "'\n";
        }
      }
      // Now add to overlap_links for legacy code, but also store extended info for mobius_links below.
      for (const auto &ol : overlap_links_ext)
        overlap_links.push_back({ol.ov_ball, ol.dst_ball});
    }
  }

  // Cross-set Möbius pairs: for each ball B with a non-self dest_ball D,
  // emit one entry per type-neighbour of B.  During solving we enforce that
  // M (recomputed from current ball positions every 50 iterations) maps the
  // source conformal vector σ̂_s to the dest conformal vector σ̂_d:
  //   e = M·σ̂_s − σ̂_d = 0
  // This is a direct 5D residual constraint on M itself, far more precise
  // than the old pairwise inversive-distance proxy.
  struct MobiusPair { Set::Ball *src; Set::Ball *dst; int src_ti; int dst_ti; };
  std::vector<MobiusPair> mobius_pairs;

  // Deduplicated list of unique src→dst links for M recomputation.
  struct MobiusLink { Set::Ball *src; Set::Ball *dst; std::vector<int> dst_ti_for_src_ti; };
  std::vector<MobiusLink> mobius_links;
  std::set<std::tuple<Set::Ball*, Set::Ball*>> mobius_link_keys;
  auto addMobiusLink = [&](Set::Ball *src, Set::Ball *dst) {
    if (src == nullptr || dst == nullptr || src == dst) return;
    auto key = std::make_tuple(src, dst);
    if (mobius_link_keys.insert(key).second)
      mobius_links.push_back({src, dst, identityTypeMap(*src, *dst)});
  };

  for (auto &set : sets)
  {
    for (auto &ball : set.balls)
    {
      if (ball.dest_ball == nullptr || ball.dest_ball == &ball) continue;
      addMobiusLink(&ball, ball.dest_ball);
    }
  }

  // For overlap_links, we need to know which parent mobius to use for the src ball.
  // So we use the extended struct if available.
  for (const auto &lk : overlap_links)
    addMobiusLink(lk.first, lk.second); // legacy, for now

  // Extended: for overlap_links_ext, store which parent mobius to use for each link.
  // We'll use this info in the mobius_pairs loop below.
  static std::vector<std::tuple<Set::Ball*, Set::Ball*, int>> mobius_parent_for_overlap;
  mobius_parent_for_overlap.clear();
  for (const auto &ol : overlap_links_ext)
    mobius_parent_for_overlap.push_back({ol.ov_ball, ol.dst_ball, ol.parent_mobius});

  // Build per-neighbour constraints from all (regular + overlap) links.
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
  // If false: keep cross-set constraints, but fix every link transform to identity.
  const bool fit_mobius_transform = true;
  if (!fit_mobius_transform)
  {
    for (auto &lk : mobius_links)
      lk.src->mobius = Set::Ball::Mobius();
  }

  std::mt19937 rng(42);
  bool warmup_done = !fit_mobius_transform;

  // Proximity anchors: snapshot initial spheres and softly keep the solve near them.
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

  // ── Step 3: iterate ──────────────────────────────────────────────────────
  // The first `warmup_iters` iterations run intra constraints only.  Once the
  // per-set configurations are reasonably converged we compute M from those
  // positions and hold it FIXED for the remainder of the solve.
  //
  // M must NOT be refreshed during the joint solve.  Doing so creates a
  // collapsing attractor: as src and dst drift toward each other the refreshed
  // M→I, which then tightens the src≈dst constraint further — ending with
  // both sets at the same position rather than a proper Möbius image of each.
  const int warmup_iters = 1500;
  const int effective_warmup_iters = fit_mobius_transform ? warmup_iters : 0;
  for (int it = 0; it < iterations + effective_warmup_iters; it++)
  {
    // After the warm-up phase, compute M once and start applying mobius pairs.
    if (fit_mobius_transform && !warmup_done && it >= warmup_iters)
    {
      for (auto &lk : mobius_links)
      {
        std::cout << "fitting Mobius transform" << std::endl;
        lk.dst_ti_for_src_ti = findBestNeighbourTypeMap(*lk.src, *lk.dst);
        computeMobiusTransform(*lk.src, *lk.dst, /*quiet=*/true, &lk.dst_ti_for_src_ti);
        std::cout << "fitted Mobius transform" << std::endl;
      }
      rebuildMobiusPairs();
      warmup_done = true;
    }

    // Shuffle both pools independently to avoid ordering bias.
    std::shuffle(intra_pairs.begin(), intra_pairs.end(), rng);
    std::shuffle(mobius_pairs.begin(), mobius_pairs.end(), rng);

    // ── (A) Within-set constraints ─────────────────────────────────────────
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
          // Rescue non-intersecting / containing pairs toward tangency.
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
    // ── (B) Möbius transform residual constraints ─────────────────────────
    //
    // For each type-neighbour ti of each src→dst link, enforce M·σ̂_s = σ̂_d
    // where M = src.mobius.M (refreshed every 50 iterations) and
    //   σ̂(C,r) = (C/r, (1−|C|²+r²)/(2r), (1+|C|²−r²)/(2r))
    //
    // Objective: f = ½‖e‖²  where  e = M·σ̂_s − σ̂_d.
    //   ∂f/∂σ̂_s = Mᵀe,  ∂f/∂σ̂_d = −e.
    //
    // Chain dσ̂/d(C,r):
    //   gC = (gs[0:3] + (gs[4]−gs[3])·C) / r
    //   gr = (−C/r²)·gs[0:3] + (r²−1+|C|²)/(2r²)·gs[3] − (r²+1+|C|²)/(2r²)·gs[4]
    //
    // GS step: δ = −(½‖e‖²/‖∇f‖²)·∇f  — zeros f in one linear step.
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

    if (!warmup_done) continue; // only apply cross-set constraints after warm-up

    for (auto &mp : mobius_pairs)
    {
      Set::Ball &src_ball = *mp.src;
      Set::Ball &dst_ball = *mp.dst;

      int si = src_ball.type_to_set[mp.src_ti];
      int di = dst_ball.type_to_set[mp.dst_ti];

      Set::Ball &sA = src_ball.parent_set->balls[si];
      Set::Ball &dA = dst_ball.parent_set->balls[di];

      double rS = clampSignedRadius(sA.radius), rD = clampSignedRadius(dA.radius);
      Eigen::Vector3d CS = sA.centre, CD = dA.centre;

      // Check if this is an overlap constraint and if so, which parent mobius to use.
      // Default: use src_ball.mobius.M
      const Set::Ball::Mobius *mobius_to_use = &src_ball.mobius;
      for (const auto &tup : mobius_parent_for_overlap) {
        if (std::get<0>(tup) == &sA && std::get<1>(tup) == &dA) {
          // Use parent0 or parent1 mobius from the overlap's parent set.
          int parent_idx = std::get<2>(tup);
          const Set *parent_set = sA.parent_set;
          const auto &ov = parent_set->overlaps;
          // Find the overlap struct that matches this sA.
          for (const auto &ovr : ov) {
            if ((parent_idx == 0 && &parent_set->balls[ovr.ball_0] == &sA) ||
                (parent_idx == 1 && &parent_set->balls[ovr.ball_1] == &sA)) {
              // Use the correct parent's mobius.
              if (parent_idx == 0)
                mobius_to_use = &parent_set->balls[ovr.ball_0].mobius;
              else
                mobius_to_use = &parent_set->balls[ovr.ball_1].mobius;
              break;
            }
          }
          break;
        }
      }

      Vec5 sigma_s = conformal_sphere(CS, rS);
      Vec5 sigma_d = conformal_sphere(CD, rD);

      // Residual: how far M·σ̂_s is from σ̂_d.
      Vec5 e = mobius_to_use->M * sigma_s - sigma_d;

      // Gradient of ½‖e‖² w.r.t. σ̂_s is Mᵀe; w.r.t. σ̂_d is −e.
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

    // ── (C) Formal proximity anchors ─────────────────────────────────────
    // Adds weighted constraints that keep each ball close to its initial
    // centre/radius while still allowing motion to satisfy topology.
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

        // Gradient descent on anchor_weight * (||C-C0||^2 + (r-r0)^2)/2.
        b.centre -= (anchor_step * w) * (b.centre - C0);
        b.radius = clampSignedRadius(b.radius - (anchor_step * w) * (b.radius - r0));
      }
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
