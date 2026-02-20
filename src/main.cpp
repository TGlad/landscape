#include <iostream>
#include "landscape.h"

int main()
{
  // Example: A simple 3-node path 0-1-2
  std::vector<int> G1 = {0, 1, 0, 0, 1, 0};
  Adj g1;
  g1.data = G1;
  std::vector<int> G2 = {0, 1, 0, 1, 0, 0};
  Adj g2;
  g2.data = G2;
  std::string hash1 = g1.get_wl_hash();
  std::string hash2 = g2.get_wl_hash();

  std::cout << "Graph 1 Hash: " << hash1 << std::endl;
  std::cout << "Graph 2 Hash: " << hash2 << std::endl;

  if (hash1 == hash2) {
      std::cout << "The graphs are likely isomorphic!" << std::endl;
  } else {
      std::cout << "The graphs are definitely different." << std::endl;
  }

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
