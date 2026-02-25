#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

static const double pi = std::acos(-1.0);

#include <iostream>
#include <vector>
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

void Landscape::calculateLeafBalls()
{
  for (auto &set: sets)
  {
    for (auto &lbi: set.leaf_ball_ids)
    {
      if (lbi[3] == -1)
        set.calculateLeafBall(lbi[0], lbi[1], lbi[2]);
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
                                   Landscape::Set::Ball &dst,
                                   bool quiet = false)
{
  src.mobius.M = Mat5::Identity();

  int m = (int)src.type_to_set.size(); // index 0 = ball itself
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
        std::cout << " ball " << i << ":";
        Set *last_par = nullptr;
        for (auto *dest_ball: ball.type->balls)
        {
          Set *par = dest_ball->parent_set;
          for (int j = 0; j<par->balls.size(); j++)
          {
            if (&par->balls[j] == dest_ball)
            {
              if (par != last_par)
                std::cout << " " << par->name << ": ";
              std::cout << j << ",";
            }
          }
          last_par = par;
        }
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

  // ── Step 2: build all constraint pairs ──────────────────────────────────
  //
  // Each entry is one scalar constraint.  We use a tagged union approach:
  // either a within-set pair (conn order), or a cross-set Gram-matching pair.
  const double damping = 1e-10;
  const double k = 0.1; // minimum extra gap for separation constraints

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

  // Cross-set Möbius pairs: for each ball B with a non-self dest_ball D,
  // emit one entry per type-neighbour of B.  During solving we enforce that
  // M (recomputed from current ball positions every 50 iterations) maps the
  // source conformal vector σ̂_s to the dest conformal vector σ̂_d:
  //   e = M·σ̂_s − σ̂_d = 0
  // This is a direct 5D residual constraint on M itself, far more precise
  // than the old pairwise inversive-distance proxy.
  struct MobiusPair { Set::Ball *src; Set::Ball *dst; int ti; };
  std::vector<MobiusPair> mobius_pairs;

  // Deduplicated list of unique src→dst links for M recomputation.
  struct MobiusLink { Set::Ball *src; Set::Ball *dst; };
  std::vector<MobiusLink> mobius_links;

  for (auto &set : sets)
  {
    for (auto &ball : set.balls)
    {
      if (ball.dest_ball == nullptr || ball.dest_ball == &ball) continue;
      mobius_links.push_back({&ball, ball.dest_ball});
      int m = (int)std::min(ball.type_to_set.size(), ball.dest_ball->type_to_set.size());
      for (int ti = 0; ti < m; ti++)
        mobius_pairs.push_back({&ball, ball.dest_ball, ti});
    }
  }

  std::mt19937 rng(42);
  bool warmup_done = false;

  // ── Step 3: iterate ──────────────────────────────────────────────────────
  // The first `warmup_iters` iterations run intra constraints only.  Once the
  // per-set configurations are reasonably converged we compute M from those
  // positions and hold it FIXED for the remainder of the solve.
  //
  // M must NOT be refreshed during the joint solve.  Doing so creates a
  // collapsing attractor: as src and dst drift toward each other the refreshed
  // M→I, which then tightens the src≈dst constraint further — ending with
  // both sets at the same position rather than a proper Möbius image of each.
  const int warmup_iters = 500;
  for (int it = 0; it < iterations + warmup_iters; it++)
  {
    // After the warm-up phase, compute M once and start applying mobius pairs.
    if (!warmup_done && it >= warmup_iters)
    {
      for (auto &lk : mobius_links)
        computeMobiusTransform(*lk.src, *lk.dst, /*quiet=*/true);
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

      // Fixed balls contribute nothing to the gradient and receive no update.
      if (bi.is_fixed) { g_dir_i = Eigen::Vector3d::Zero(); g_dist_i = 0; g_curv_i = 0; }
      if (bj.is_fixed) { g_dir_j = Eigen::Vector3d::Zero(); g_dist_j = 0; g_curv_j = 0; }

      double g2 = g_dir_i.squaredNorm() + g_dist_i*g_dist_i + g_curv_i*g_curv_i
                + g_dir_j.squaredNorm() + g_dist_j*g_dist_j + g_curv_j*g_curv_j;
      double step = -error / (g2 + damping);
      bi.dir       += step * g_dir_i;  bi.dir.normalize();
      bi.dist      += step * g_dist_i;
      if (bi.curvature != 0.0) bi.curvature = std::max(1e-6, bi.curvature + step * g_curv_i);
      bj.dir       += step * g_dir_j;  bj.dir.normalize();
      bj.dist      += step * g_dist_j;
      if (bj.curvature != 0.0) bj.curvature = std::max(1e-6, bj.curvature + step * g_curv_j);
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
    // Chain dσ̂/d(C,r) (Jacobian rows: I/r | −Cᵀ/r | Cᵀ/r) then
    // d(C,r)/d(dir,dist,curvature):
    //   gC = (gs[0:3] + (gs[4]−gs[3])·C) / r
    //   gr = (−C/r²)·gs[0:3] + (r²−1+|C|²)/(2r²)·gs[3] − (r²+1+|C|²)/(2r²)·gs[4]
    //   g_dir  = (dist+r)·(gC − (gC·d̂)·d̂)
    //   g_dist = gC·d̂
    //   g_curv = (gC·d̂ + gr)·(−1/κ²)
    //
    // GS step: δ = −(½‖e‖²/‖∇f‖²)·∇f  — zeros f in one linear step.
    auto sigma_to_ball_grad = [](const Set::Ball &b, double r,
                                 const Eigen::Vector3d &C, const Vec5 &gs)
        -> std::tuple<Eigen::Vector3d, double, double>
    {
      double C2 = C.squaredNorm();
      Eigen::Vector3d gC = (gs.head<3>() + (gs(4) - gs(3)) * C) / r;
      double gr = (-C / (r*r)).dot(gs.head<3>())
                + (r*r - 1.0 + C2) / (2.0*r*r) * gs(3)
                - (r*r + 1.0 + C2) / (2.0*r*r) * gs(4);
      Eigen::Vector3d g_dir  = (b.dist + r) * (gC - gC.dot(b.dir) * b.dir);
      double          g_dist = gC.dot(b.dir);
      double          g_curv = (gC.dot(b.dir) + gr) * (-1.0 / (b.curvature * b.curvature));
      return {g_dir, g_dist, g_curv};
    };

    if (!warmup_done) continue; // only apply cross-set constraints after warm-up

    for (auto &mp : mobius_pairs)
    {
      Set::Ball &src_ball = *mp.src;
      Set::Ball &dst_ball = *mp.dst;

      int si = src_ball.type_to_set[mp.ti];
      int di = dst_ball.type_to_set[mp.ti];

      Set::Ball &sA = src_ball.parent_set->balls[si];
      Set::Ball &dA = dst_ball.parent_set->balls[di];

      if (sA.curvature == 0.0 || dA.curvature == 0.0) continue;

      double rS = 1.0/sA.curvature, rD = 1.0/dA.curvature;
      Eigen::Vector3d CS = sA.dir*(sA.dist+rS), CD = dA.dir*(dA.dist+rD);

      Vec5 sigma_s = conformal_sphere(CS, rS);
      Vec5 sigma_d = conformal_sphere(CD, rD);

      // Residual: how far M·σ̂_s is from σ̂_d.
      Vec5 e = src_ball.mobius.M * sigma_s - sigma_d;

      // Gradient of ½‖e‖² w.r.t. σ̂_s is Mᵀe; w.r.t. σ̂_d is −e.
      Vec5 g_sigma_s = src_ball.mobius.M.transpose() * e;
      Vec5 g_sigma_d = -e;

      auto [gDirS, gDistS, gCurvS] = sigma_to_ball_grad(sA, rS, CS, g_sigma_s);
      auto [gDirD, gDistD, gCurvD] = sigma_to_ball_grad(dA, rD, CD, g_sigma_d);

      if (sA.is_fixed) { gDirS = Eigen::Vector3d::Zero(); gDistS = 0.0; gCurvS = 0.0; }
      if (dA.is_fixed) { gDirD = Eigen::Vector3d::Zero(); gDistD = 0.0; gCurvD = 0.0; }

      double g2 = gDirS.squaredNorm() + gDistS*gDistS + gCurvS*gCurvS
                + gDirD.squaredNorm() + gDistD*gDistD + gCurvD*gCurvD;

      double step = -e.squaredNorm() / (2.0 * (g2 + damping));

      sA.dir += step * gDirS;  sA.dir.normalize();
      sA.dist += step * gDistS;
      if (sA.curvature != 0.0) sA.curvature = std::max(1e-6, sA.curvature + step * gCurvS);
      dA.dir += step * gDirD;  dA.dir.normalize();
      dA.dist += step * gDistD;
      if (dA.curvature != 0.0) dA.curvature = std::max(1e-6, dA.curvature + step * gCurvD);
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