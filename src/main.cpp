#include <iostream>
#include "landscape.h"

int main()
{
  // OK let's make some sets manually to start off with.
  Landscape land;
  land.sets.push_back(Landscape::Set("tree", 6));
  Landscape::Set &tree = land.sets.back();
 /*     struct Ball // supports oriented spheres and oriented planes
    {
      Eigen::Vector3d dir;
      double dist;
      double curvature;
      int sets_id;  // indexes sets list
      int types_id; // indexes types list
      int balls_id; // indexes Type's sphere list
    };*/
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    tree.conn(0, i) = 3;
    tree.conn(5, i) = 3;
    // around the meridian
    tree.conn(last_i, i) = 3;
    last_i = i;
  }
  // now come up with some approximate locations
  tree.balls[0].dir = Eigen::Vector3d(0,0,1);
  tree.balls[1].dir = Eigen::Vector3d(1,0,0);
  tree.balls[2].dir = Eigen::Vector3d(0,1,0);
  tree.balls[3].dir = Eigen::Vector3d(-1,0,0);
  tree.balls[4].dir = Eigen::Vector3d(0,-1,0);
  tree.balls[5].dir = Eigen::Vector3d(0,0,-1);
  for (int i = 0; i<6; i++)
  {
    tree.balls[i].dist = 1.0;
    tree.balls[i].curvature = 1.0;
  }
  // type and dest_ball need to get figured out later by accumulating the types as we go

  tree.applyConnectivity();
  land.addSetToTypes(tree);

  return 0;
}
