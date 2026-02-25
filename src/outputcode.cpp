#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

static const double pi = std::acos(-1.0);

#include <iostream>
#include <vector>
#include <Eigen/Dense>

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

  // Find the flat index of a Ball pointer into BALLS[]
  auto flatIndex = [&](const Set::Ball *ball) -> int {
    for (int si = 0; si < num_sets; si++)
      for (int bi = 0; bi < (int)sets[si].balls.size(); bi++)
        if (&sets[si].balls[bi] == ball) return offsets[si] + bi;
    return -1;
  };

  // Find which set index a Ball pointer belongs to (-1 if null)
  auto setIndex = [&](const Set::Ball *ball) -> int {
    if (!ball) return -1;
    for (int si = 0; si < num_sets; si++)
      for (int bi = 0; bi < (int)sets[si].balls.size(); bi++)
        if (&sets[si].balls[bi] == ball) return si;
    return -1;
  };

  out << std::fixed << std::setprecision(7);

  // ── Ball struct ──────────────────────────────────────────────────────────
  out << "struct Ball {\n"
      << "    vec3  dir;\n"
      << "    float dist;\n"
      << "    float curvature;\n"
      << "    int   dest_set;   // flat index into SETS[];  -1 = reflexive\n"
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
      int dest_s = setIndex(b.dest_ball);
      int dest_b = flatIndex(b.dest_ball);
      bool last = (si == num_sets - 1 && bi == n - 1);
      out << "    Ball(vec3("
          << b.dir.x() << ", " << b.dir.y() << ", " << b.dir.z() << "), "
          << b.dist << ", "
          << b.curvature << ", "
          << dest_s << ", "
          << dest_b << ")"
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
            << "-1, -1)"  // leaf balls don't recurse
            << (flat < total_leaf_balls - 1 ? "," : "") << "\n";
      }
    }
    out << ");\n\n";
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
        << "vec3 applyMobius(int mi, vec3 v, inout float scale) {\n"
        << "    if (MOBIUS_SIM[mi])\n"
        << "    {\n"
        << "        scale *= MOBIUS_S[mi];\n"
        << "        return MOBIUS_T[mi] + MOBIUS_S[mi] * (MOBIUS_R[mi] * v);\n"
        << "    }\n"
        << "    vec3 w = v - MOBIUS_C[mi];\n"
        << "    scale *= MOBIUS_S[mi] / dot(w, w);\n"
        << "    return MOBIUS_T[mi] + (MOBIUS_S[mi] / dot(w, w)) * (MOBIUS_R[mi] * w);\n"
        << "}\n";
  }

  std::cout << "Wrote " << filename << "\n";
}
