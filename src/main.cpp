#include <iostream>
#include "landscape.h"

int main()
{
  // OK let's make some sets manually to start off with.
  Landscape land;

  land.sets.push_back(Landscape::Set("tree-tree", 6));
  Landscape::Set &tree_tree = land.sets.back();
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    tree_tree.conn(0, i) = 3;
    tree_tree.conn(5, i) = 3;
    // around the meridian
    tree_tree.conn(last_i, i) = 3;
    last_i = i;
  }
  // now come up with some approximate locations
  tree_tree.balls[0].dir = Eigen::Vector3d(0,0,1);
  tree_tree.balls[1].dir = Eigen::Vector3d(1,0,0);
  tree_tree.balls[2].dir = Eigen::Vector3d(0,1,0);
  tree_tree.balls[3].dir = Eigen::Vector3d(-1,0,0);
  tree_tree.balls[4].dir = Eigen::Vector3d(0,-1,0);
  tree_tree.balls[5].dir = Eigen::Vector3d(0,0,-1);
  for (int i = 0; i<6; i++)
  {
    tree_tree.balls[i].dist = 1.0;
    tree_tree.balls[i].curvature = 1.0;
  }

  tree_tree.applyConnectivity();
  land.addSetToTypes(tree_tree);


  land.sets.push_back(Landscape::Set("shell-tree", 6));
  Landscape::Set &shell_tree = land.sets.back();
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    shell_tree.conn(0, i) = 3;
    shell_tree.conn(last_i, i) = i==1 ? 2 : 3;
    last_i = i;
  }
  shell_tree.conn(1, 3) = 5;
  shell_tree.conn(1, 5) = 3;
  shell_tree.conn(2, 5) = 3;
  shell_tree.conn(3, 5) = 2;
  
  // now come up with some approximate locations
  shell_tree.balls[0].dir = Eigen::Vector3d(0,0,1);
  shell_tree.balls[1].dir = Eigen::Vector3d(0.5,-1,0).normalized();
  shell_tree.balls[2].dir = Eigen::Vector3d(0.5, 1,0).normalized();
  shell_tree.balls[3].dir = Eigen::Vector3d(-1,0,0);
  shell_tree.balls[4].dir = Eigen::Vector3d(0.3,0,0.3).normalized();
  shell_tree.balls[5].dir = Eigen::Vector3d(0,0,-1);
  for (int i = 0; i<6; i++)
  {
    shell_tree.balls[i].dist = 1.0;
    shell_tree.balls[i].curvature = 1.0;
  }
  shell_tree.balls[4].dist = 0.3;
  shell_tree.balls[4].curvature = 8.0;

  shell_tree.applyConnectivity();
  land.addSetToTypes(shell_tree);
  // only once all sets have been added can we try and match up the destination names and numbers

  land.matchUpDestinationBalls();
  return 0;
}
