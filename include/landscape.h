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
    Set(const std::string &name, int num_balls);
    std::string name;
    Eigen::Vector4d colour {Eigen::Vector4d(1,1,1,1)};
    struct Ball // supports oriented spheres and oriented planes
    {
      Eigen::Vector3d dir;
      double dist;
      double curvature;
      std::string dest_set;
      int dest_ball_id {-1}; // -1 will pick the first in the set that works
      bool is_fixed {false};

      void initSphere(const Eigen::Vector3d &p, double rad, bool fix = false);
      void initPlane(const Eigen::Vector3d &normal, double d, bool fix = false);

      // auto-set
      Set *parent_set; 
      Type *type;
      Ball *dest_ball {nullptr};
      std::vector<int> type_to_set; // type_to_set[canonical_type_idx] = set ball index
      std::vector<int> set_to_type; // set_to_type[set_ball_idx]  = canonical_type_idx (-1 if not neighbour) 
      struct Mobius
      {
        Eigen::Matrix<double,5,5> M = Eigen::Matrix<double,5,5>::Identity();
        // decomposed version of M:
        Eigen::Vector3d C = Eigen::Vector3d::Zero();
        Eigen::Vector3d T = Eigen::Vector3d::Zero();
        double          s = 1.0;
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
        bool is_similarity = true;

        Eigen::Vector3d transformPoint(const Eigen::Vector3d &p) const;
        Eigen::Vector3d transformDecomposed(const Eigen::Vector3d &p) const;
        std::pair<Eigen::Vector3d,double> transformSphere( const Eigen::Vector3d &c, double r) const;
      } mobius;
    };
    std::vector<Ball> balls;
    std::vector<Ball> leaf_balls; // used to represent set at leaf

    Adj conn; // connectivity. -1=kissing, 0 is disconnected
    bool leaf_union; // union if true, else intersection
    bool render_volume_only {false};

    void addLeafBall(int i, int j, int k, int l);
    void addLeafBall(int i, int j, int k); // smallest: center in plane of 3 ball centers
    void findOrthogonalSphere(int I, int J, int K, int L);

    // internal stuff
    std::vector<Eigen::Vector4i> leaf_ball_ids;
    bool verifyConnectivity(double tol = 1e-4) const;
    void calculateLeafBall(int i, int j, int k, int l);
    void calculateLeafBall(int i, int j, int k); // smallest: center in plane of 3 ball centers

    // Copy/move constructors must re-point parent_set to *this, not to the source.
    Set(const Set &o);
    Set(Set &&o);
    Set &operator=(const Set &o);
    Set &operator=(Set &&o);
  };
  std::deque<Set> sets;

  struct Type
  {
    Type(const Adj &set_connectivity, int ball_id, std::vector<int> &out_is);
    Adj conn;
    std::vector<Set::Ball *> balls; // sets id, balls id
  };
  std::deque<Type> types;

  void printConnectivity(bool show_valid_destinations = false); // called after all addSetToTypes() are called
  void applyConnectivity(int iterations = 4000);
  void verifyConnectivity();  
  void calculateLeaves();
  void addSetsToTypes();  
  void addSetToTypes(Set &set);
  void matchUpDestinationBalls();
  void outputCode(const std::string &filename = "landscape.glsl") const;
};