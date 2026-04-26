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
      std::vector<int> location; // sequence of ids to pass through before swapping to dest_set 
      int dest_ball_id; // -1 choses the first that works topologically
      double mobility {1.0}; // 0=fixed, 1=standard update weight

      void initSphere(const Eigen::Vector3d &p, double rad, double mobility_value = 1.0);
      void initPlane(const Eigen::Vector3d &normal, double d, double mobility_value = 1.0);

      inline double radius() const { return 1.0/curvature; }
      inline Eigen::Vector3d centre() const { return dir * (dist + radius());}
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
    std::vector<int> leaf_ball_set;
    struct Overlap // supports oriented spheres and oriented planes
    {
      Overlap(int ball0, int ball1, const std::string &destset, int dest_ball_id) : ball_0(ball0), ball_1(ball1), dest_set(destset), dest_ball_1_id(dest_ball_id) {}
      int ball_0, ball_1;
      std::string dest_set;
      int dest_ball_1_id; // dest_ball_0_id is balls[ball_0].dest_id
    };
    std::vector<Overlap> overlaps;


    Adj conn; // connectivity. -1=kissing, 0 is disconnected
    bool leaf_union; // union if true, else intersection
    bool render_volume_only {false};

    void addLeafBall(int i, int j, int k, int l);
    void addLeafBall(int i, int j, int k, double scale = 1.0); // for n=3 (bald) points
    void addLeafBalls(std::vector<int> leaf_ids){  leaf_ball_set = leaf_ids; }
    void findOrthogonalSphere(int I, int J, int K, int L);

    // internal stuff
    std::vector<Eigen::Vector4i> leaf_ball_ids;
    std::vector<double> leaf_ball_scales;
    bool verifyConnectivity(double tol = 1e-4) const;
    void calculateLeafBall(int i, int j, int k, int l);
    void calculateLeafBall(int i, int j, int k, double scale); // smallest: center in plane of 3 ball centers
    void calculateLeafBalls(); // uses leaf_ball_ids

    // Copy/move constructors must re-point parent_set to *this, not to the source.
    Set(const Set &o);
    Set(Set &&o);
    Set &operator=(const Set &o);
    Set &operator=(Set &&o);
  };
  std::deque<Set> sets;
  const Set &set(const std::string &set_name)
  {
    for (const auto &set: sets)
    {
      if (set.name == set_name)
        return set;
    }
    Set set("dummy", 0);
    return set;
  }

  struct Type
  {
    Type(const Adj &set_connectivity, int ball_id, std::vector<int> &out_is);
    Adj conn;
    std::vector<Set::Ball *> balls; // sets id, balls id
  };
  std::deque<Type> types;

  void printConnectivity(bool show_valid_destinations = false); // called after all addSetToTypes() are called
  void applyConnectivity(int iterations = 6000);
  void verifyConnectivity();  
  void calculateLeaves();
  void addSetsToTypes();  
  void addSetToTypes(Set &set);
  void matchUpDestinationBalls();
  void outputCode(const std::string &filename = "landscape.glsl") const;
};
