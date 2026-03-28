auto cube_sphere = [&]()
{
  land.sets.push_back(Landscape::Set("cube_sphere", 5));
  Landscape::Set &set = land.sets.back();

  // between spheres
  set.conn(0,1) = 3; 

  // between planes
  set.conn(2,4) = 2; // 45 degrees
  set.conn(3,4) = 4; // 45 degrees
  set.conn(2,3) = 3; 

  // spheres to planes
  set.conn(0,2) = 0; // or 6?
  set.conn(0,3) = 2; 
  set.conn(0,4) = 2; 
  set.conn(1,2) = 2; 
  set.conn(1,3) = 2; 
  set.conn(1,4) = 6; // since between balls would be 3

  set.balls[0].initSphere(Eigen::Vector3d(0.8, 0.0, 0.0), 0.3);
  set.balls[1].initSphere(Eigen::Vector3d(1.0, 1.0, 1.0).normalized(), 2.0/3.0);
  set.balls[2].initPlane(Eigen::Vector3d(1,-1,0),  0.0);
  set.balls[3].initPlane(Eigen::Vector3d(0,1,-1),  0.0);
  set.balls[4].initPlane(Eigen::Vector3d(0,0, 1),  0.0);

//  set.addLeafBall(0,2,3,4);
  set.addLeafBall(1,2,3,4);
  set.render_volume_only = true;   

  set.balls[0].dest_set = "cube_sphere2";
  set.balls[0].dest_ball_id = 0;
};

auto cube_sphere2 = [&]()
{
  land.sets.push_back(Landscape::Set("cube_sphere2", 5));
  Landscape::Set &set = land.sets.back();

  // between spheres
  set.conn(0,1) = 3; 

  // between planes
  set.conn(2,4) = 2; // 45 degrees
  set.conn(3,4) = 4; // 45 degrees
  set.conn(2,3) = 3; 

  // spheres to planes
  set.conn(0,2) = 0; // or 6?
  set.conn(0,3) = 2; 
  set.conn(0,4) = 2; 
  set.conn(1,2) = 2; 
  set.conn(1,3) = 2; 
  set.conn(1,4) = 6; // since between balls would be 3

  set.balls[0].initSphere(Eigen::Vector3d(0.8, 0.0, 0.0), 0.3);
  set.balls[1].initSphere(Eigen::Vector3d(1.0, 1.0, 1.0).normalized(), 2.0/3.0);
  set.balls[2].initPlane(Eigen::Vector3d(1,-1,0),  0.0);
  set.balls[3].initPlane(Eigen::Vector3d(0,1,-1),  0.0);
  set.balls[4].initPlane(Eigen::Vector3d(0,0, 1),  0.0);

//  set.addLeafBall(0,2,3,4);
  set.addLeafBall(1,2,3,4);
  set.render_volume_only = true;   

  set.balls[0].dest_set = "cube_tree";
  set.balls[0].dest_ball_id = 0;
};

auto cube_tree = [&]()
{
  land.sets.push_back(Landscape::Set("cube_tree", 5));
  Landscape::Set &set = land.sets.back();

  // between spheres
  set.conn(0,1) = 3; 

  // between planes
  set.conn(2,4) = 2; // 45 degrees
  set.conn(3,4) = 4; // 45 degrees
  set.conn(2,3) = 3; 

  // spheres to planes
  set.conn(0,2) = 0; // or 6?
  set.conn(0,3) = 2; 
  set.conn(0,4) = 2; 
  set.conn(1,2) = 2; 
  set.conn(1,3) = 2; 
  set.conn(1,4) = 6; // since between balls would be 3

  set.balls[0].initSphere(Eigen::Vector3d(1.5, 0.0, 0.0), 0.6);
  set.balls[1].initSphere(Eigen::Vector3d(1.0, 1.0, 1.0).normalized(), 2.0/3.0);
  set.balls[2].initPlane(Eigen::Vector3d(1,-1,0),  0.0);
  set.balls[3].initPlane(Eigen::Vector3d(0,1,-1),  0.0);
  set.balls[4].initPlane(Eigen::Vector3d(0,0, 1),  0.0);

//  set.addLeafBall(0,2,3,4);
  set.addLeafBall(1,2,3,4);
  set.render_volume_only = true;   

//  set.balls[0].dest_set = "ball1";
//  set.balls[0].dest_ball_id = 0;
};

auto ball = [&]()
{
  land.sets.push_back(Landscape::Set("ball", 6));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.8,0.6,0.5,1);
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
  // now come up with some approximate locations
  constexpr double noise = 0.5;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };
  
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0,1)), 1.0);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d(1,0,0)), 1.0);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(0,1,0)), 1.0);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(-1,0,0)), 1.0);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,-1,0)), 1.0);
  set.balls[5].initSphere(with_noise(Eigen::Vector3d(0,0,-1)), 1.0);

//  set.balls[0].location.push_back(0); // go through ball 5 first
//  set.balls[0].location.push_back(5); // go through ball 5 first
//  set.balls[0].dest_set = "shell-shell";
//  set.balls[0].dest_ball_id = 0;
  set.addLeafBall(0,1,2);
  set.addLeafBall(0,2,3);
  set.addLeafBall(0,3,4);
  set.addLeafBall(0,4,1);
  set.addLeafBall(5,1,2);
  set.addLeafBall(5,2,3);
  set.addLeafBall(5,3,4);
  set.addLeafBall(5,4,1);
//  set.addLeafBall(0,1,2,4);
//  set.addLeafBall(2,3,4,5);

  set.leaf_union = true;
  set.render_volume_only = true; 
};

auto clusterTree2 = [&]()
{
  land.sets.push_back(Landscape::Set("cluster-tree2", 10));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.1,0.8,0.1,1);
  // like octahedron of order 3
  // but additional sphere...
  int last_i = 4;
  for (int i = 1; i<=4; i++)
  {
    set.conn(0, i) = 3;
    set.conn(5, i) = 3;
    // around the meridian
    set.conn(last_i, i) = 3;
    last_i = i;
  }
  set.conn(6, 0) = 3;
  set.conn(6, 1) = 3;
  set.conn(6, 2) = 3;
  set.conn(7, 0) = 3;
  set.conn(7, 2) = 3;
  set.conn(7, 3) = 3;
  set.conn(8, 0) = 3;
  set.conn(8, 3) = 3;
  set.conn(8, 4) = 3;
  set.conn(9, 0) = 3;
  set.conn(9, 4) = 3;
  set.conn(9, 1) = 3;
  // now come up with some approximate locations
  constexpr double noise = 0.2;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };
  
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0,1)), 1.0);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d(1,0,0)), 1.0);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(0,1,0)), 1.0);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(-1,0,0)), 1.0);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,-1,0)), 1.0);
  set.balls[5].initSphere(with_noise(Eigen::Vector3d(0,0,-1)), 1.0);
  set.balls[6].initSphere(with_noise(Eigen::Vector3d(1,1,1)), 1.0);
  set.balls[7].initSphere(with_noise(Eigen::Vector3d(-1,1,1)), 0.8);
  set.balls[8].initSphere(with_noise(Eigen::Vector3d(-1,-1,1)), 0.5);
  set.balls[9].initSphere(with_noise(Eigen::Vector3d(1,-1,1)), 1.3);

//  set.balls[0].location.push_back(0); // go through ball 5 first
//  set.balls[0].location.push_back(5); // go through ball 5 first
//  set.balls[0].dest_set = "shell-shell";
//  set.balls[0].dest_ball_id = 0;
  set.addLeafBall(0,1,2);
  set.addLeafBall(0,2,3);
  set.addLeafBall(0,3,4);
  set.addLeafBall(0,4,1);
  set.addLeafBall(5,1,2);
  set.addLeafBall(5,2,3);
  set.addLeafBall(5,3,4);
  set.addLeafBall(5,4,1);
  set.addLeafBall(0,1,2,6); // tetrahedron at end
  set.addLeafBall(0,2,3,7);
  set.addLeafBall(0,3,4,8);
  set.addLeafBall(0,4,1,9);

  set.leaf_union = true;
  set.render_volume_only = true; 
};