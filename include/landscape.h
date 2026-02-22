#pragma once

#include <string>
#include <deque>
#include <numeric>
#include "adj.h"
#include <Eigen/Dense>

struct Landscape
{
  struct Type; // forward declaration so Set::Ball::type resolves to Landscape::Type
  struct Set
  {
    Set(const std::string &name, int num_balls) : name(name) 
    { 
      balls.resize(num_balls); 
      for (int i = 0; i<num_balls; i++)
      {
        balls[i].parent_set = this;
      }
      conn.resize(num_balls); // defaults to disconnected
    }
    std::string name;
    struct Ball // supports oriented spheres and oriented planes
    {
      Eigen::Vector3d dir;
      double dist;
      double curvature;
      std::string dest_set;
      int dest_ball_id {-1}; // -1 will pick the first in the set that works

      // auto-set
      Set *parent_set; 
      Type *type;
      Ball *dest_ball {nullptr};
      std::vector<int> type_to_set; // type_to_set[canonical_type_idx] = set ball index
      std::vector<int> set_to_type; // set_to_type[set_ball_idx]  = canonical_type_idx (-1 if not neighbour) 

      struct Mobius
      {
        Eigen::Vector3d center {Eigen::Vector3d(0,0,0)};
        Eigen::Vector3d translation {Eigen::Vector3d(0,0,0)};
        Eigen::Matrix3d rotation {Eigen::Matrix3d::Identity()};
        double scale {1.0};
        bool flip {false};

        Eigen::Vector3d transformMobius(const Eigen::Vector3d &p)
        {
          Eigen::Vector3d v = p - center;
          double dist_sq = flip ? v.squaredNorm() : 1.0;
          Eigen::Vector3d v_inv = v / dist_sq;
          return (rotation * (scale * v_inv)) + translation;
        }
      } mobius;
    };
    std::vector<Ball> balls;
    std::vector<Ball> leaf_balls; // used to represent set at leaf
    Adj conn; // connectivity. -1=kissing, 0 is disconnected
    bool leaf_union; // union if true, else intersection
    void applyConnectivity();
    bool verifyConnectivity(double tol = 1e-4) const;
    void addLeafBall(int i, int j, int k, int l);
    void addLeafBall(int i, int j, int k); // smallest: center in plane of 3 ball centers
  };
  std::deque<Set> sets;

  struct Type
  {
    Type(const Adj &set_connectivity, int ball_id, std::vector<int> &out_is)
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
    Adj conn;
    std::vector<Set::Ball *> balls; // sets id, balls id
  };
  std::deque<Type> types;

  void addSetToTypes(Set &set);

  void matchUpDestinationBalls();
  void outputCode(const std::string &filename = "landscape.glsl") const;
};