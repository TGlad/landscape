#pragma once

#include <string>
#include <deque>
#include "adj.h"
#include <Eigen/Dense>

struct Landscape
{
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
      Set *parent_set; 
      struct Type *type;      
      int dest_ball; // indexes Type's sphere list
    };
    std::vector<Ball> balls;
    Adj conn; // connectivity. -1=kissing, 0 is disconnected
    void applyConnectivity();
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

  void addSetToTypes(Set &set)
  {
    for (int i = 0; i<(int)set.balls.size(); i++)
    { 
      Type new_type(set.conn, i);
      bool found = false;
      for (int j = 0; j<types.size(); j++)
      {
        if (types[j].conn == new_type.conn)
        {
          // new_type already exists, so add it in.
          types[j].balls.push_back(&set.balls[i]);
          found = true;
          break;
        }
      }
      if (!found)
      {
        types.push_back(new_type);
        types.back().balls.push_back(&set.balls[i]);
      }
    }
  }
};