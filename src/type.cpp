#include "landscape.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>

Landscape::Type::Type(const Adj &set_connectivity, int ball_id, std::vector<int> &out_is)
{
  int n = set_connectivity.size();
  std::vector<int> is;
  is.push_back(ball_id);
  for (int i = 0; i<n; i++)
  {
    if (set_connectivity(i, ball_id) > 0 && i!=ball_id)
      is.push_back(i);
  }
  int m = (int)is.size();
  conn.resize(m);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      conn(i,j) = std::max(0, set_connectivity(is[i], is[j])); // kissing points ignored

  // Canonicalize: permute indices 1..m-1 to get the lex-smallest lower triangle,
  // so that conn can be compared directly regardless of original neighbour ordering.
  if (m > 2)
  {
    // Build the lower-triangle key for a given permutation of slots 1..m-1.
    // Slot 0 stays fixed; perm[k] is the conn-row used for slot k+1.
    auto make_key = [&](const std::vector<int> &perm) {
      std::vector<int> key;
      key.reserve(m * (m + 1) / 2);
      for (int i = 1; i < m; i++)
      {
        key.push_back(conn(perm[i-1], 0)); // connection to node 0
        for (int j = 1; j < i; j++)
          key.push_back(conn(perm[i-1], perm[j-1]));
      }
      return key;
    };

    std::vector<int> perm(m - 1);
    std::iota(perm.begin(), perm.end(), 1); // [1, 2, ..., m-1]
    std::vector<int> best_perm = perm;
    std::vector<int> best_key  = make_key(perm);

    while (std::next_permutation(perm.begin(), perm.end()))
    {
      auto key = make_key(perm);
      if (key < best_key) { best_key = key; best_perm = perm; }
    }

    // Rebuild conn using the canonical permutation.
    Adj canonical;
    canonical.resize(m);
    for (int i = 0; i < m; i++)
    {
      int ri = (i == 0) ? 0 : best_perm[i-1];
      for (int j = 0; j <= i; j++)
      {
        int rj = (j == 0) ? 0 : best_perm[j-1];
        canonical(i, j) = conn(ri, rj);
      }
    }
    conn = std::move(canonical);

    // Write out the final canonical is[] ordering.
    out_is.resize(m);
    out_is[0] = is[0];
    for (int i = 1; i < m; i++)
      out_is[i] = is[best_perm[i-1]];
  }
  else
  {
    out_is = is; // m <= 2: no permutation needed
  }
}
