#include <iostream>
#include <iomanip>
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
    // adjacent around the meridian contact at pi/3
    tree_tree.conn(last_i, i) = 3;
    last_i = i;
  }
  // opposing equatorial balls are disjoint (explicitly separate)
  tree_tree.conn(1, 3) = 0;
  tree_tree.conn(2, 4) = 0;
  // now come up with some approximate locations
  double eq = 0.3; // equatorial distance and radius
  tree_tree.balls[0].initSphere(Eigen::Vector3d(0,0,1), 1.0);
  tree_tree.balls[1].initSphere(Eigen::Vector3d(eq,0,0), eq);
  tree_tree.balls[2].initSphere(Eigen::Vector3d(0,eq,0), eq);
  tree_tree.balls[3].initSphere(Eigen::Vector3d(-eq,0,0), eq);
  tree_tree.balls[4].initSphere(Eigen::Vector3d(0,-eq,0), eq);
  tree_tree.balls[5].initSphere(Eigen::Vector3d(0,0,-1), 1.0);
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
  shell_tree.balls[0].initSphere(Eigen::Vector3d(0,0,1), 1.0);
  shell_tree.balls[1].initSphere(Eigen::Vector3d(0.5,-1,0), 1.0);
  shell_tree.balls[2].initSphere(Eigen::Vector3d(0.5, 1,0), 1.0);
  shell_tree.balls[3].initSphere(Eigen::Vector3d(-1,0,0), 1.0);
  shell_tree.balls[4].initSphere(Eigen::Vector3d(0.3,0,0.3), 0.125);
  shell_tree.balls[5].initSphere(Eigen::Vector3d(0,0,-1), 1.0);
  land.addSetToTypes(shell_tree);

  land.sets.push_back(Landscape::Set("shell-shell", 6));
  Landscape::Set &shell_shell = land.sets.back();
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    shell_shell.conn(0, i) = 3;
    shell_shell.conn(5, i) = 3;
    // around the meridian
    shell_shell.conn(last_i, i) = 3;
    last_i = i;
  }
  // now come up with some approximate locations
  shell_shell.balls[0].initSphere(Eigen::Vector3d(0,0,0.35), 0.25);
  shell_shell.balls[1].initSphere(Eigen::Vector3d(1,0,0), 1.0);
  shell_shell.balls[2].initSphere(Eigen::Vector3d(0,1,0), 1.0);
  shell_shell.balls[3].initSphere(Eigen::Vector3d(-1,0,0), 1.0);
  shell_shell.balls[4].initSphere(Eigen::Vector3d(0,-1,0), 1.0);
  shell_shell.balls[5].initSphere(Eigen::Vector3d(0,0,-0.35), 0.25);

  shell_shell.balls[0].dest_set = "tree-tree";
  shell_shell.balls[0].dest_ball_id = 0;
  land.addSetToTypes(shell_shell);


  land.sets.push_back(Landscape::Set("shell-sponge", 6));
  Landscape::Set &shell_sponge = land.sets.back();
  shell_sponge = shell_shell; // shortcut
  shell_sponge.name = "shell-sponge"; // restore name after copy
  shell_sponge.conn(0,5) = 4; // make it a sponge
  shell_sponge.balls[0].initSphere(Eigen::Vector3d(0,0,0.251), 0.25);
  shell_sponge.balls[5].initSphere(Eigen::Vector3d(0,0,-0.251), 0.25);
  shell_sponge.balls[0].dest_set = "";
  land.addSetToTypes(shell_sponge);

  land.sets.push_back(Landscape::Set("cube_tree-tree", 4));
  Landscape::Set &cube_tree_tree = land.sets.back();
 // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  last_i = 3;
  for (int i = 1; i<=3; i++)
  {
    cube_tree_tree.conn(0, i) = 6;
    cube_tree_tree.conn(last_i, i) = 2;
    last_i = i;
  }
  cube_tree_tree.balls[0].initSphere(Eigen::Vector3d(1,1,1), 1.2);
  cube_tree_tree.balls[1].initPlane(Eigen::Vector3d(1,0,0),  0.0, true);
  cube_tree_tree.balls[2].initPlane(Eigen::Vector3d(0,1,0),  0.0, true);
  cube_tree_tree.balls[3].initPlane(Eigen::Vector3d(0,0,1),  0.0, true);
  land.addSetToTypes(cube_tree_tree);

  // ── Global joint solve ───────────────────────────────────────────────────
  // Solves all sets simultaneously: per-set connectivity constraints AND
  // the inversive-distance matching constraints imposed by dest_set links.
  land.applyConnectivity();
  land.printConnectivity();

  tree_tree.verifyConnectivity();
  shell_tree.verifyConnectivity();
  shell_shell.verifyConnectivity();
  shell_sponge.verifyConnectivity();
  cube_tree_tree.verifyConnectivity();

  tree_tree.addLeafBall(0,1,2,3);
  tree_tree.addLeafBall(2,3,4,5);
  tree_tree.leaf_union = true;

  shell_shell.addLeafBall(0,1,2,4);
  shell_shell.addLeafBall(2,3,4,5);
  shell_shell.leaf_union = false;

  shell_sponge.addLeafBall(0,1,2,4);
  shell_sponge.addLeafBall(2,3,4,5);
  shell_sponge.leaf_union = false;

//  cube_tree_tree.findOrthogonalSphere(0,1,2,3);
 // cube_tree_tree.addLeafBall(0,1,2,3); // no leaf ball as the spheres don't overlap in 3s. We should check this automatically
  land.matchUpDestinationBalls();

  land.outputCode();
  return 0;
}
