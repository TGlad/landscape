#include <iostream>
#include <iomanip>
#include "landscape.h"

int main()
{
  Landscape land;

  #include "simple.hpp" // define some simple sets

  // decide which shapes to add:
  // trees:
  voidTree(); clusterTree(); treeTree(); shellTree();
  // sponges
  voidSponge(); clusterSponge(); spongeSponge(); shellSponge();
  // shells
  shellShell();

  // misc
  cubeSpongeSponge();

  #include "items.hpp" // define some simple sets
  cube_sphere();
  cube_sphere2();
  cube_tree();
  ball();

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
  return 0;
}
