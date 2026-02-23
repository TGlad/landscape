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
  tree_tree.balls[0].dir = Eigen::Vector3d(0,0,1);
  tree_tree.balls[1].dir = Eigen::Vector3d(1,0,0);
  tree_tree.balls[2].dir = Eigen::Vector3d(0,1,0);
  tree_tree.balls[3].dir = Eigen::Vector3d(-1,0,0);
  tree_tree.balls[4].dir = Eigen::Vector3d(0,-1,0);
  tree_tree.balls[5].dir = Eigen::Vector3d(0,0,-1);
  
  tree_tree.balls[0].dist = 1.0;
  tree_tree.balls[0].curvature = 1.0;
  tree_tree.balls[5].dist = 1.0;
  tree_tree.balls[5].curvature = 1.0;
  for (int i = 1; i<5; i++)
  {
    tree_tree.balls[i].dist = 0.3;
    tree_tree.balls[i].curvature = 3.0;
  }
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
  shell_shell.balls[0].dir = Eigen::Vector3d(0,0,1);
  shell_shell.balls[1].dir = Eigen::Vector3d(1,0,0);
  shell_shell.balls[2].dir = Eigen::Vector3d(0,1,0);
  shell_shell.balls[3].dir = Eigen::Vector3d(-1,0,0);
  shell_shell.balls[4].dir = Eigen::Vector3d(0,-1,0);
  shell_shell.balls[5].dir = Eigen::Vector3d(0,0,-1);
  for (int i = 0; i<6; i++)
  {
    shell_shell.balls[i].dist = 1.0;
    shell_shell.balls[i].curvature = 1.0;
  }
  shell_shell.balls[0].dist = 0.1;
  shell_shell.balls[5].dist = 0.1;
  shell_shell.balls[0].curvature = 4.0;
  shell_shell.balls[5].curvature = 4.0;
  shell_shell.balls[0].dest_set = "tree-tree";
  shell_shell.balls[0].dest_ball_id = 0;
  land.addSetToTypes(shell_shell);


  land.sets.push_back(Landscape::Set("shell-sponge", 6));
  Landscape::Set &shell_sponge = land.sets.back();
  shell_sponge = shell_shell; // shortcut
  shell_sponge.name = "shell-sponge"; // restore name after copy
  shell_sponge.conn(0,5) = 4; // make it a sponge
  shell_sponge.balls[0].dist = 0.01;
  shell_sponge.balls[5].dist = 0.01;
  shell_sponge.balls[0].dest_set = "";
  land.addSetToTypes(shell_sponge);

  // ── Global joint solve ───────────────────────────────────────────────────
  // Solves all sets simultaneously: per-set connectivity constraints AND
  // the inversive-distance matching constraints imposed by dest_set links.
  land.applyConnectivity();

  tree_tree.verifyConnectivity();
  shell_tree.verifyConnectivity();
  shell_shell.verifyConnectivity();
  shell_sponge.verifyConnectivity();

  tree_tree.addLeafBall(0,1,2,3);
  tree_tree.addLeafBall(2,3,4,5);
  tree_tree.leaf_union = true;

  shell_shell.addLeafBall(0,1,2,4);
  shell_shell.addLeafBall(2,3,4,5);
  shell_shell.leaf_union = false;

  shell_sponge.addLeafBall(0,1,2,4);
  shell_sponge.addLeafBall(2,3,4,5);
  shell_sponge.leaf_union = false;

  land.matchUpDestinationBalls();

  land.outputCode();
  return 0;
}
