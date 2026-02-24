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
    std::string name;
    struct Ball // supports oriented spheres and oriented planes
    {
      Eigen::Vector3d dir;
      double dist;
      double curvature;
      std::string dest_set;
      int dest_ball_id {-1}; // -1 will pick the first in the set that works
      bool is_fixed {false};

      void initSphere(const Eigen::Vector3d &p, double rad)
      {
        dir = p.normalized();
        dist = p.norm() - rad;
        curvature = 1.0/rad;
      }
      void initPlane(const Eigen::Vector3d &normal, double d)
      {
        dir = normal.normalized();
        dist = d;
        curvature = 0.0;
      }

      // auto-set
      Set *parent_set; 
      Type *type;
      Ball *dest_ball {nullptr};
      std::vector<int> type_to_set; // type_to_set[canonical_type_idx] = set ball index
      std::vector<int> set_to_type; // set_to_type[set_ball_idx]  = canonical_type_idx (-1 if not neighbour) 

      struct Mobius
      {
        // O(4,1) matrix kept for solver use (Gauss-Seidel, etc.)
        // Metric η = diag(1,1,1,1,−1).
        // Point p encoded as (p, (1−|p|²)/2, (1+|p|²)/2); decode Y → Y.xyz/(Y[3]+Y[4]).
        Eigen::Matrix<double,5,5> M = Eigen::Matrix<double,5,5>::Identity();

        // GLSL-friendly decomposition, populated from M by decomposeMobius().
        //
        //   Similarity   (is_similarity = true,  C unused):
        //     f(v) = T + s * R * v
        //
        //   Inversion    (is_similarity = false):
        //     f(v) = T + s * R * (v − C) / |v − C|²
        //     C = pre-image of ∞,  T = image of ∞
        Eigen::Vector3d C = Eigen::Vector3d::Zero();
        Eigen::Vector3d T = Eigen::Vector3d::Zero();
        double          s = 1.0;
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        bool is_similarity = true;

        Eigen::Vector3d transformPoint(const Eigen::Vector3d &p) const
        {
          double r2 = p.squaredNorm();
          Eigen::Matrix<double,5,1> P;
          P << p.x(), p.y(), p.z(), (1.0-r2)/2.0, (1.0+r2)/2.0;
          auto Y = M * P;
          double w = Y(3) + Y(4);
          return (std::abs(w) > 1e-15) ? Eigen::Vector3d(Y(0)/w, Y(1)/w, Y(2)/w)
                                       : Eigen::Vector3d::Zero();
        }

        // Apply the T,C,s,R decomposition (GLSL path) to a point.
        //   Similarity (is_similarity=true):  f(v) = T + s*R*v
        //   Inversion  (is_similarity=false): f(v) = T + s*R*(v-C)/|v-C|²
        Eigen::Vector3d transformDecomposed(const Eigen::Vector3d &p) const
        {
          if (is_similarity)
            return T + s * (R * p);
          Eigen::Vector3d w = p - C;
          double d2 = w.squaredNorm();
          if (d2 < 1e-30) return Eigen::Vector3d(1e15, 0.0, 0.0);
          return T + (s / d2) * (R * w);
        }

        // Returns (centre, radius) of the transformed sphere.
        std::pair<Eigen::Vector3d,double> transformSphere(
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
      } mobius;
    };
    std::vector<Ball> balls;
    std::vector<Ball> leaf_balls; // used to represent set at leaf
    Adj conn; // connectivity. -1=kissing, 0 is disconnected
    bool leaf_union; // union if true, else intersection

    bool verifyConnectivity(double tol = 1e-4) const;
    void addLeafBall(int i, int j, int k, int l);
    void addLeafBall(int i, int j, int k); // smallest: center in plane of 3 ball centers

    Set(const std::string &name, int num_balls) : name(name) 
    { 
      balls.resize(num_balls); 
      for (int i = 0; i<num_balls; i++)
        balls[i].parent_set = this;
      conn.resize(num_balls); // defaults to disconnected
    }
    // Copy/move constructors must re-point parent_set to *this, not to the source.
    Set(const Set &o)
      : name(o.name), balls(o.balls), leaf_balls(o.leaf_balls),
        conn(o.conn), leaf_union(o.leaf_union)
    {
      for (auto &b : balls)      b.parent_set = this;
      for (auto &b : leaf_balls) b.parent_set = this;
    }
    Set(Set &&o)
      : name(std::move(o.name)), balls(std::move(o.balls)),
        leaf_balls(std::move(o.leaf_balls)),
        conn(std::move(o.conn)), leaf_union(o.leaf_union)
    {
      for (auto &b : balls)      b.parent_set = this;
      for (auto &b : leaf_balls) b.parent_set = this;
    }
    Set &operator=(const Set &o)
    {
      if (this == &o) return *this;
      name = o.name; balls = o.balls; leaf_balls = o.leaf_balls;
      conn = o.conn; leaf_union = o.leaf_union;
      for (auto &b : balls)      b.parent_set = this;
      for (auto &b : leaf_balls) b.parent_set = this;
      return *this;
    }
    Set &operator=(Set &&o)
    {
      if (this == &o) return *this;
      name = std::move(o.name); balls = std::move(o.balls);
      leaf_balls = std::move(o.leaf_balls);
      conn = std::move(o.conn); leaf_union = o.leaf_union;
      for (auto &b : balls)      b.parent_set = this;
      for (auto &b : leaf_balls) b.parent_set = this;
      return *this;
    }
  };
  std::deque<Set> sets;

  void printConnectivity(bool show_valid_destinations = false); // called after all addSetToTypes() are called
  // Joint Gauss-Seidel over all sets simultaneously.
  // Enforces each set's own pairwise connectivity constraints AND, for every
  // ball that has a dest_set link, the inversive-distance matching constraints
  // that ensure a Möbius transform exists between the two neighbourhoods.
  // Call this once after all sets (and their dest_set/dest_ball_id fields)
  // are configured, instead of calling Set::applyConnectivity() individually.
  void applyConnectivity(int iterations = 4000);

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