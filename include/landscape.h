#pragma once

#include <string>
#include <deque>
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
      Ball *dest_ball; 
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
    Type(const Adj &set_connectivity, int ball_id)
    {
      int n = set_connectivity.size();
      std::vector<int> is;
      for (int i = 0; i<n; i++)
      {
        if (set_connectivity(i, ball_id) > 0 || i==ball_id)
          is.push_back(i);
      }
      conn.resize((int)is.size());
      
      for (int i = 0; i<(int)is.size(); i++)
      {
        int I = is[i];
        for (int j = 0; j<(int)is.size(); j++)
        {
          int J = is[j];
          conn(i,j) = std::max(0, set_connectivity(I, J)); // kissing points ignored
        }
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