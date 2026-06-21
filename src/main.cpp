#include <iostream>
#include <iomanip>
#include "landscape.h"

int main()
{
  Landscape land;

  #include "simple.hpp" // define some simple sets

//  clusterTree();
//  shellShell();
  icosahedron();

  #include "items.hpp" // define some simple sets
//  ball();
//  clusterTree2();
//  bush1();
//  bush2();
//  bush3();
//  ballToLandscape();
  tree_test();
  hill_test();
 // hill_testb();

  // ── Global joint solve ───────────────────────────────────────────────────
  // Solves all sets simultaneously: per-set connectivity constraints AND
  // the inversive-distance matching constraints imposed by dest_set links.
  land.generateOverlapLayouts();

  // users can modify overlap layouts here if they want something different

  land.addSetsToTypes();
  land.applyConnectivity();
  land.verifyConnectivity();
//  land.printConnectivity(true);
  land.calculateLeaves();
  land.matchUpDestinationBalls();
  land.outputCode();

/*
  std::vector<int> ids = {0,1,2,3,4,5,6,10};
  auto &s1 = land.set("icosahedron");
  auto &s2 = land.set("tree-test");
  auto &s3 = land.set("hill-testb");
  auto &s4 = land.set("tree-hill-test");

  std::cout << "ball6: " << s4.balls[6].radius << " and ball10: " << s2.balls[10].radius << " should equal one of: " << s2.balls[6].radius << ", " << s2.balls[7].radius << ", " << s2.balls[8].radius << ", " << s2.balls[9].radius << ", " << s2.balls[10].radius << ", " << std::endl;
  std::cout << "ball3: " << s4.balls[3].radius << " and ball4: " << s4.balls[4].radius<< " should equal one of: " << s3.balls[6].radius << ", " << s3.balls[7].radius << ", " << s3.balls[8].radius << ", " << s3.balls[9].radius << ", " << s3.balls[10].radius << std::endl;
*//*  for (auto &id: ids)
  {
    std::cout << id << " centres: " << s1.balls[id].centre.transpose() << " : " << s2.balls[id].centre.transpose() << " : " << s3.balls[id].centre.transpose() << " : " << s4.balls[id].centre.transpose() << ", radii: " << s1.balls[id].radius << " : " << s2.balls[id].radius << " : " << s3.balls[id].radius << " : " << s4.balls[id].radius << std::endl;
  }*/
  return 0;
}
