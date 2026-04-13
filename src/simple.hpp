auto voidTree = [&]()
{
  land.sets.push_back(Landscape::Set("void-tree", 7));
  Landscape::Set &set = land.sets.back();
  for (int i = 0; i<7; i++)
  {
    if (i != 3)
      set.conn(i, 3) = -1;
  }
  set.balls[0].initSphere(Eigen::Vector3d(0,0,1), 0.8);
  set.balls[1].initSphere(Eigen::Vector3d(1,0,0), 0.8);
  set.balls[2].initSphere(Eigen::Vector3d(1,0,0), 0.8);
  set.balls[3].initSphere(Eigen::Vector3d(0,0,0), 0.2);
  set.balls[4].initSphere(Eigen::Vector3d(0,1,0), 0.8);
  set.balls[5].initSphere(Eigen::Vector3d(0,-1,0), 0.8);
  set.balls[6].initSphere(Eigen::Vector3d(0,0,-1), 0.8);
};

auto voidSponge = [&]()
{
  land.sets.push_back(Landscape::Set("void-sponge", 5));
  Landscape::Set &set = land.sets.back();
  for (int i = 0; i<5; i++)
    for (int j = i+1; j<5; j++)
      set.conn(i, j) = 4;
  set.balls[0].initSphere(Eigen::Vector3d(0.6,0,-0.3), 0.8);
  set.balls[1].initSphere(Eigen::Vector3d(-0.3,-0.6,-0.3), 0.8);
  set.balls[2].initSphere(Eigen::Vector3d(-0.3,0.6,-0.3), 0.8);
  set.balls[3].initSphere(Eigen::Vector3d(0,0,0.6), 0.8);
  set.balls[4].initSphere(Eigen::Vector3d(0,0,0), 0.2);
};

auto clusterTree = [&]()
{
  land.sets.push_back(Landscape::Set("cluster-tree", 5));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.25,0.7,0.2,1);
  // double tetrahedron
  int last_i = 3;
  for (int i = 1; i<=3; i++)
  {
    set.conn(0, i) = 3;
    set.conn(4, i) = 3;
    set.conn(last_i, i) = 3;
    last_i = i;
  }
  // now come up with some approximate locations
  double eq = 0.3; // equatorial distance and radius
  set.balls[0].initSphere(Eigen::Vector3d(0,0.3, 1), 0.8);
  set.balls[1].initSphere(Eigen::Vector3d( 1,-1, 0), 0.8);
  set.balls[2].initSphere(Eigen::Vector3d(-1,-1, 0), 0.8);
  set.balls[3].initSphere(Eigen::Vector3d(0,0,0), 0.2);
  set.balls[4].initSphere(Eigen::Vector3d(0,0.3,-1), 0.8);
  set.addLeafBall(0,1,2,3);
  set.addLeafBall(1,2,3,4);
  set.leaf_union = true;
};

auto clusterSponge = [&]()
{
  land.sets.push_back(Landscape::Set("cluster-sponge", 5));
  Landscape::Set &set = land.sets.back();
  for (int i = 0; i<5; i++)
    for (int j = i+1; j<5; j++)
      set.conn(i, j) = 3;
  for (int i = 1; i<5; i++)
    set.conn(0,i) = 4;
  set.balls[0].initSphere(Eigen::Vector3d(0,0,0.6), 0.8);
  set.balls[1].initSphere(Eigen::Vector3d(0.6,0,-0.3), 0.8);
  set.balls[2].initSphere(Eigen::Vector3d(-0.3,-0.6,-0.3), 0.8);
  set.balls[3].initSphere(Eigen::Vector3d(-0.3,0.6,-0.3), 0.8);
  set.balls[4].initSphere(Eigen::Vector3d(0,0,0), 0.2);
  set.addLeafBall(1,2,3,4);
  set.render_volume_only = true; // optional?
};

auto treeTree = [&]()
{
  land.sets.push_back(Landscape::Set("tree-tree", 6));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0,1,0,1);
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    set.conn(0, i) = 3;
    set.conn(5, i) = 3;
    // adjacent around the meridian contact at pi/3
    set.conn(last_i, i) = 3;
    last_i = i;
  }
  // now come up with some approximate locations
  double eq = 0.25; // equatorial distance and radius
  set.balls[0].initSphere(Eigen::Vector3d(0, 0, 1), 1.0);
  set.balls[1].initSphere(Eigen::Vector3d(eq,0, 0), eq);
  set.balls[2].initSphere(Eigen::Vector3d(0,eq, 0), eq);
  set.balls[3].initSphere(Eigen::Vector3d(-eq,0,0), eq);
  set.balls[4].initSphere(Eigen::Vector3d(0,-eq,0), eq);
  set.balls[5].initSphere(Eigen::Vector3d(0, 0,-1), 1.0);
  set.addLeafBall(0,1,2,3);
  set.addLeafBall(2,3,4,5);
  set.leaf_union = true;
  set.render_volume_only = true; 
};

auto shellTree = [&]()
{
  land.sets.push_back(Landscape::Set("shell-tree", 6));
  Landscape::Set &set = land.sets.back();
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    set.conn(0, i) = 3;
    set.conn(last_i, i) = i==1 ? 2 : 3;
    last_i = i;
  }
  set.conn(1, 3) = 5;
  set.conn(1, 5) = 3;
  set.conn(2, 5) = 3;
  set.conn(3, 5) = 2;
  // now come up with some approximate locations
  set.balls[0].initSphere(Eigen::Vector3d(0,0,1), 1.0);
  set.balls[1].initSphere(Eigen::Vector3d(0.5,-1,0), 1.0);
  set.balls[2].initSphere(Eigen::Vector3d(0.5, 1,0), 1.0);
  set.balls[3].initSphere(Eigen::Vector3d(-1,0,0), 1.0);
  set.balls[4].initSphere(Eigen::Vector3d(0.3,0,0.3), 0.125);
  set.balls[5].initSphere(Eigen::Vector3d(0,0,-1), 1.0);
};

auto shellShell = [&]()
{
  land.sets.push_back(Landscape::Set("shell-shell", 6));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.7,0.5,0.35,1);
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    set.conn(0, i) = 3;
    set.conn(5, i) = 3;
    // around the meridian
    set.conn(last_i, i) = 3;
    last_i = i;
  }

  constexpr double noise = 0.0;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };    

  // now come up with some approximate locations
  double r = std::sqrt(2.0/3.0);
  double m = 0.5;
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0,0.5)), 0.5);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d(1,0,0)), r, 1.0);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(0,1,0)), r, m);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(-1,0,0)), r, 1.0);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,-1,0)), r, m);
  set.balls[5].initSphere(with_noise(Eigen::Vector3d(0,0,-0.5)), 0.5);

//    set.balls[5].dest_set = "cluster-tree2";
//    set.balls[5].dest_ball_id = 5;
  set.addLeafBall(0,1,2,4);
  set.addLeafBall(2,3,4,5);
  set.leaf_union = false;
//    set.addLeafBalls({0,1,2,3,4,5});
  set.render_volume_only = true; 
};

auto spongeSponge = [&]()
{
  land.sets.push_back(Landscape::Set("sponge-sponge", 6));
  Landscape::Set &set = land.sets.back();
  // octahedron based
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    set.conn(0, i) = 3;
    set.conn(5, i) = 3;
    // around the meridian
    set.conn(last_i, i) = 3;
    last_i = i;
  }
  set.conn(0,5) = 4; // make it a sponge
  // now come up with some approximate locations
  set.balls[0].initSphere(Eigen::Vector3d(-0.5,0,1), 1.1);
  set.balls[1].initSphere(Eigen::Vector3d(0,-0.05,0), 0.06);
  set.balls[2].initSphere(Eigen::Vector3d(0, 0.05,0), 0.06);
  set.balls[3].initSphere(Eigen::Vector3d(0.27,-0.25,0), 0.3);
  set.balls[4].initSphere(Eigen::Vector3d(0.27,0.25,0), 0.3);
  set.balls[5].initSphere(Eigen::Vector3d(-0.5,0,-1), 1.1);
  set.addLeafBall(0,1,2,4);
  set.addLeafBall(2,3,4,5);
  set.leaf_union = true;
};

auto shellSponge = [&]()
{
  land.sets.push_back(Landscape::Set("shell-sponge", 6));
  Landscape::Set &set = land.sets.back();
  // octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    set.conn(0, i) = 3;
    set.conn(5, i) = 3;
    // around the meridian
    set.conn(last_i, i) = 3;
    last_i = i;
  }
  set.conn(0,5) = 4; // make it a sponge
  // now come up with some approximate locations
  set.balls[0].initSphere(Eigen::Vector3d(0,0,0.251), 0.25);
  set.balls[1].initSphere(Eigen::Vector3d(1,0,0), 1.0);
  set.balls[2].initSphere(Eigen::Vector3d(0,1,0), 1.0);
  set.balls[3].initSphere(Eigen::Vector3d(-1,0,0), 1.0);
  set.balls[4].initSphere(Eigen::Vector3d(0,-1,0), 1.0);
  set.balls[5].initSphere(Eigen::Vector3d(0,0,-0.251), 0.25);
  set.addLeafBall(0,1,2,4);
  set.addLeafBall(2,3,4,5);
  set.leaf_union = false;
  set.balls[2].dest_set = "sponge-sponge";
  set.balls[2].dest_ball_id = 3;
};

auto cubeSpongeSponge = [&]()
{
  land.sets.push_back(Landscape::Set("cube_sponge-sponge", 4));
  Landscape::Set &set = land.sets.back();
// octahedron of order 3, so 0=top, 1,2,3,4 is mid and 5 is base
  int last_i = 3;
  for (int i = 1; i<=3; i++)
  {
    set.conn(0, i) = 5;
    set.conn(last_i, i) = 2;
    last_i = i;
  }
  set.balls[0].initSphere(Eigen::Vector3d(1,1,1), 1.2);
  set.balls[1].initPlane(Eigen::Vector3d(1,0,0),  0.0, true);
  set.balls[2].initPlane(Eigen::Vector3d(0,1,0),  0.0, true);
  set.balls[3].initPlane(Eigen::Vector3d(0,0,1),  0.0, true);
};

auto icosahedron = [&]()
{
  land.sets.push_back(Landscape::Set("icosahedron", 12));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.7,0.6,0.55,1);
  const double pi = 3.14159265;

  double h = std::sqrt(1.0 / 5.0);
  double r = std::sqrt(1.0 - h); // 0.743
  set.balls[0].initSphere(Eigen::Vector3d(0,0,1), r);
  for (int i = 0; i<5; i++)
  {
    float ang1 = (double)i * 2.0*pi/5.0;
    float ang2 = ang1 + pi/5.0;

    set.conn(0,1+i) = 2; // top fan
    set.balls[1+i].initSphere(h*Eigen::Vector3d(2.0*std::cos(ang1),2.0*std::sin(ang1), 1), r, 0.0); // top ring
    set.conn(1+i, 1 + (i+1)%5) = 2; // around top ring
    set.conn(1+i, 6+i) = 2; // zig
    set.conn(6+i, 1+ (i+1)%5) = 2; // zag
    set.balls[6+i].initSphere(h*Eigen::Vector3d(2.0*std::cos(ang2),2.0*std::sin(ang2),-1), r); // bottom ring
    set.conn(6+i, 6 + (i+1)%5) = 2; // around bottom ring
    set.conn(6+i,11) = 2; // bottom fan
  }
  set.balls[11].initSphere(Eigen::Vector3d(0,0,-1), r);

  set.balls[0].dest_set = "tree-test";
  set.balls[0].dest_ball_id = 0;
  set.balls[0].overlap_sets.resize(2);
  set.balls[0].overlap_ids.resize(2, -1);
  set.balls[0].overlap_sets[1] = "tree-hill-test";
  set.balls[0].overlap_ids[1] = 1; // <-- index set_ball_id should not be an overlap destination
 // set.balls[1].location.push_back(1);
  set.balls[1].dest_set = "hill-test";
  set.balls[1].dest_ball_id = 1;

  set.ball_pair[0][1].set_set = "tree-hill-test";
  set.ball_pair[]


  set.addLeafBalls({0,1,2,4,5,6,7,8,9,10,11});
  set.leaf_union = false;
  set.render_volume_only = true; 
};
