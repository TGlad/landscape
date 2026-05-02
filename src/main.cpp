#include <iostream>
#include <iomanip>
#include "landscape.h"

int main()
{
  Landscape land;

  #include "simple.hpp" // define some simple sets

  // decide which shapes to add:
  // trees:
  //voidTree(); clusterTree(); treeTree(); shellTree();
  // sponges
//  voidSponge(); clusterSponge(); spongeSponge(); shellSponge();
  // shells
//  shellShell();
  icosahedron();
  icos2();

  // misc
//  cubeSpongeSponge();

  #include "items.hpp" // define some simple sets
//  cube_sphere();
//  cube_sphere2();
//  cube_tree();
//  ball();
//  clusterTree2();
//  bush1();
//  bush2();
//  bush3();
//  ballToLandscape();
  tree_test();
  hill_test();
  tree_hill_test();
//  hill_icos_test();

  // ── Global joint solve ───────────────────────────────────────────────────
  // Solves all sets simultaneously: per-set connectivity constraints AND
  // the inversive-distance matching constraints imposed by dest_set links.
  land.addSetsToTypes();
  land.applyConnectivity();
  land.verifyConnectivity();
  land.printConnectivity(true);
  land.calculateLeaves();
  land.matchUpDestinationBalls();
  land.outputCode();


/*  std::vector<int> ids = {0,1,2,3,4,5,6,10};
  auto &s1 = land.set("icosahedron");
  auto &s2 = land.set("tree-test");
  auto &s3 = land.set("hill-test");
  auto &s4 = land.set("tree-hill-test");
  for (auto &id: ids)
  {
    std::cout << id << " centres: " << s1.balls[id].centre.transpose() << " : " << s2.balls[id].centre.transpose() << " : " << s3.balls[id].centre.transpose() << " : " << s4.balls[id].centre.transpose() << ", radii: " << s1.balls[id].radius << " : " << s2.balls[id].radius << " : " << s3.balls[id].radius << " : " << s4.balls[id].radius << std::endl;
  }*/
  return 0;
}
