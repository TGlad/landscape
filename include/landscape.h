#pragma once

#include <string>
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
      conn.resize((num_balls*(num_balls+1))/2);
    }
    std::string name;
    struct Ball // supports oriented spheres and oriented planes
    {
      Eigen::Vector3d dir;
      double dist;
      double curvature;
      Set *parent_set;  // indexes sets list
      int type;      // indexes types list
      int dest_ball; // indexes Type's sphere list
    };
    std::vector<Ball> balls;
    std::vector<int> conn; // connectivity
    int &connection(int i, int j) { return i >= j ? conn[(i * (i + 1) / 2) + j] : conn[(j * (j + 1) / 2) + i]; }  
    void applyConnectivity();
  };
  std::vector<Set> sets;

  struct Type
  {
    Type(std::vector<int> set_connectivity, int ball_id)
    {

    }
    std::vector<int> conn;
    std::vector<Eigen::Vector2i> balls; // sets id, balls id
  };
  std::vector<Type> types;

  // do we tell type to pull in a new type from a set?
  // or do we tell set to add in a type?
  void addSetToTypes(Set &set)
  {
    for (int i = 0; i<(int)set.balls.size(); i++)
    { 
      Type new_type(set.connection, i);
      for (int j = 0; j<types.size(); j++)
      {
        if (types[i].conn == new_type.conn)
        {
          
        }
      }
    }
  }
};