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
  // now come up with some approximate locations
  constexpr double noise = 0.0;//0.5;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };
  
  double mobility = 0.5;
  double r = std::sqrt(2.0/3.0);
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0,1)), 1.0, mobility);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d(1,0,0)), r, 0.0);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(0,1,0)), r, 0.0);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(-1,0,0)), r, 0.0);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,-1,0)), r, 0.0);
  set.balls[5].initSphere(with_noise(Eigen::Vector3d(0,0,-1)), 1.0, mobility);

  set.balls[0].dest_set = "shell-shell";
  set.balls[0].dest_ball_id = 0;
/*  double scale = 0.5;
  set.addLeafBall(0,1,2, scale);
  set.addLeafBall(0,2,3, scale);
  set.addLeafBall(0,3,4, scale);
  set.addLeafBall(0,4,1, scale);
  set.addLeafBall(5,1,2, scale);
  set.addLeafBall(5,2,3, scale);
  set.addLeafBall(5,3,4, scale);
  set.addLeafBall(5,4,1, scale);*/

  set.addLeafBalls({0,1,2,3,4,5});


  set.leaf_union = true;
  set.render_volume_only = true; 
};

auto clusterTree2 = [&]()
{
  land.sets.push_back(Landscape::Set("cluster-tree2", 10));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.6,0.5,0.3,1);
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
  set.balls[6].initSphere(with_noise(Eigen::Vector3d(1,1,1)), 1.2);
  set.balls[7].initSphere(with_noise(Eigen::Vector3d(-1,1,1)), 1.0);
  set.balls[8].initSphere(with_noise(Eigen::Vector3d(-1,-1,1)), 0.5);
  set.balls[9].initSphere(with_noise(Eigen::Vector3d(1,-1,1)), 1.3);

  set.balls[6].dest_set = "cluster-tree";
  set.balls[6].dest_ball_id = 0;
  set.balls[7].dest_set = "bush1";
  set.balls[7].dest_ball_id = 0;
  set.balls[8].dest_set = "bush2";
  set.balls[8].dest_ball_id = 0;
  set.balls[9].dest_set = "bush3";
  set.balls[9].dest_ball_id = 0;
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

auto bush1 = [&]()
{
  land.sets.push_back(Landscape::Set("bush1", 5));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.25,0.6,0.2,1);
  // double tetrahedron
  int last_i = 3;
  for (int i = 1; i<=3; i++)
  {
    set.conn(0, i) = 3;
    set.conn(4, i) = 3;
    set.conn(last_i, i) = 3;
    last_i = i;
  }

  constexpr double noise = 0.2;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };

  // now come up with some approximate locations
  double eq = 0.2; // equatorial distance and radius
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0.3, 1)), 0.8);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d( 1,-1, 0)), 0.8);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(-1,-1, 0)), 0.8);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(0,0,0)), 0.2);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,0.3,-1)), 0.8);
  set.addLeafBall(0,1,2,3);
  set.addLeafBall(1,2,3,4);
  set.leaf_union = true;
};

auto bush2 = [&]()
{
  land.sets.push_back(Landscape::Set("bush2", 5));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.35,0.7,0.2,1);
  // double tetrahedron
  int last_i = 3;
  for (int i = 1; i<=3; i++)
  {
    set.conn(0, i) = 3;
    set.conn(4, i) = 3;
    set.conn(last_i, i) = 3;
    last_i = i;
  }

  constexpr double noise = 0.2;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };

  // now come up with some approximate locations
  double eq = 0.2; // equatorial distance and radius
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0.3, 1)), 0.8);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d( 1,-1, 0)), 0.8);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(-1,-1, 0)), 0.8);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(0,0,0)), 0.2);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,0.3,-1)), 0.8);
  set.addLeafBall(0,1,2,3);
  set.addLeafBall(1,2,3,4);
  set.leaf_union = true;
};

auto bush3 = [&]()
{
  land.sets.push_back(Landscape::Set("bush3", 5));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.3,0.65,0.35,1);
  // double tetrahedron
  int last_i = 3;
  for (int i = 1; i<=3; i++)
  {
    set.conn(0, i) = 3;
    set.conn(4, i) = 3;
    set.conn(last_i, i) = 3;
    last_i = i;
  }

  constexpr double noise = 0.2;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };

  // now come up with some approximate locations
  double eq = 0.2; // equatorial distance and radius
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0.3, 1)), 0.8);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d( 1,-1, 0)), 0.8);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(-1,-1, 0)), 0.8);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(0,0,0)), 0.2);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,0.3,-1)), 0.8);
  set.addLeafBall(0,1,2,3);
  set.addLeafBall(1,2,3,4);
  set.leaf_union = true;
};

auto ballToLandscape = [&]()
{
  // ball has a symmetric pyramid structure, which can switch directly to
  // tree-like protrusions or shell-like indentations, but not to saddle-like
  // landscapes, so we need a transition tile
  land.sets.push_back(Landscape::Set("ball-to-landscape", 7));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.7,0.6,0.5,1);
  int last_i = 4;
  // bottom pyramid
  for (int i = 1; i<=4; i++)
  {
    set.conn(0, i) = 3;
    // around the meridian
    set.conn(last_i, i) = 3;
    last_i = i;
  }
  set.conn(5,1) = 3;
  set.conn(5,2) = 2;
  set.conn(5,3) = 3;
  set.conn(6,3) = 3;
  set.conn(6,4) = 2;
  set.conn(6,1) = 3;
  set.conn(5,6) = 2;

  // now come up with some approximate locations
  constexpr double noise = 0.0;
  auto with_noise = [&](const Eigen::Vector3d &v)
  {
    return v + noise * Eigen::Vector3d::Random();
  };
  
  double r = std::sqrt(2.0/3.0);
  set.balls[0].initSphere(with_noise(Eigen::Vector3d(0,0,1)), 1.0);
  set.balls[1].initSphere(with_noise(Eigen::Vector3d(1,0,0)), r, 0.0);
  set.balls[2].initSphere(with_noise(Eigen::Vector3d(0,1,0)), r, 0.0);
  set.balls[3].initSphere(with_noise(Eigen::Vector3d(-1,0,0)), r, 0.0);
  set.balls[4].initSphere(with_noise(Eigen::Vector3d(0,-1,0)), r, 0.0);
  set.balls[5].initSphere(with_noise(Eigen::Vector3d(0,0.6,-1.0)), 1.0);
  set.balls[6].initSphere(with_noise(Eigen::Vector3d(0,-0.6,-1.0)), 1.0);

//  set.balls[0].dest_set = "shell-shell";
//  set.balls[0].dest_ball_id = 0;
/*  double scale = 0.9;
  set.addLeafBall(0,1,2, scale);
  set.addLeafBall(0,2,3, scale);
  set.addLeafBall(0,3,4, scale);
  set.addLeafBall(0,4,1, scale);

  set.addLeafBall(5,1,2, scale);
  set.addLeafBall(5,2,3, scale);
  set.addLeafBall(6,3,4, scale);
  set.addLeafBall(6,4,1, scale);
  set.addLeafBall(5,6,1, scale);
  set.addLeafBall(5,6,3, scale);*/

  set.addLeafBalls({0,1,2,3,4,5,6});

  set.leaf_union = true;
  set.render_volume_only = true; 
};
static double mid = 1.2; // 1 is normal
static double ridge = std::sqrt(1.5);
//static double ll = 0.9;


auto tree_test = [&]()
{
  land.sets.push_back(Landscape::Set("tree-test", 12));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.25,0.7,0.2,1);
  const double pi = 3.14159265;

  double h = std::sqrt(1.0 / 5.0);
  double r = std::sqrt(1.0 - h); // 0.743
  double peak = 1.0; // 1 is normal
  double offset = std::sqrt(1.0);
  double rmid = std::sqrt(mid);
  set.balls[0].initSphere(Eigen::Vector3d(0,0,1), r, 0.0);
  for (int i = 0; i<5; i++)
  {
    float ang1 = (double)i * 2.0*pi/5.0;
    float ang2 = ang1 + pi/5.0;

    set.conn(0,1+i) = 2; // top fan
    if ((1+i)==2)
    {
      double H = std::sqrt(4.0*h*h + (ll-h)*(ll-h));
      double R = std::sqrt(H*H - r*r);
      set.balls[1+i].initSphere((ll/std::sqrt(5.0))*Eigen::Vector3d(2.0*std::cos(ang1),2.0*std::sin(ang1), 1), R, 0.0); // top ring
    }
    else
      set.balls[1+i].initSphere(h*Eigen::Vector3d(2.0*std::cos(ang1),2.0*std::sin(ang1), 1), r, 0.0); // top ring
    set.conn(1+i, 1 + (i+1)%5) = 2; // around top ring
    set.conn(1+i, 6+i) = 2; // zig
    set.conn(6+i, 1+ (i+1)%5) = 2; // zag
    double scale = rmid;
    if (i == 1 || i== 4)
      scale *= ridge;
    else
      scale /= ridge;
    if(i==2 || i==3)
      scale *= offset;
    else
      scale /= offset;
    double moveability = 1.0;//(6+i)==6 || (6+i)==10 ? 0.0 : 1.0;
    set.balls[6+i].initSphere(scale*h*Eigen::Vector3d(2.0*std::cos(ang2),2.0*std::sin(ang2), -scale), r, moveability); // bottom ring
    set.conn(6+i, 6 + (i+1)%5) = 2; // around bottom ring
    set.conn(6+i,11) = 2; // bottom fan
  }
  set.balls[11].initSphere(Eigen::Vector3d(0,0,-peak), offset*offset*ridge*ridge*peak*r);

  set.addLeafBalls({0,1,2,4,5,6,7,8,9,10,11});
  set.leaf_union = false;
  set.render_volume_only = true; 
  set.dest_set = "hill-test";
  set.levels = 1;
};
static double bulge = 1.3;

auto hill_test = [&]()
{
  land.sets.push_back(Landscape::Set("hill-test", 12));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.85,0.45,0.25,1);
  const double pi = 3.14159265;

  double h = std::sqrt(1.0 / 5.0);
  double r = std::sqrt(1.0 - h); // 0.743

  std::cout << "LL: " << ll << std::endl;

  set.balls[0].initSphere(Eigen::Vector3d(0,0,1), r, 0.0);
  for (int i = 0; i<5; i++)
  {
    float ang1 = ((double)i - 0.15) * 2.0*pi/5.0;
    float ang2 = ang1 + pi/5.0;

    set.conn(0,1+i) = 2; // top fan
    double sc = (1+i)==3 || (1+i)==4 ? bulge : 1.0;
    double moveability = (1+i)==1 || (1+i)==2 || (i+i)==5 ? 0.0 : 1.0;
    if ((1+i)==2)
    {
      double H = std::sqrt(4.0*h*h + (ll-h)*(ll-h));
      double R = std::sqrt(H*H - r*r);
      set.balls[1+i].initSphere((ll/std::sqrt(5.0))*Eigen::Vector3d(2.0*std::cos(ang1),2.0*std::sin(ang1), 1), R, 0.0); // top ring
    }
    else
      set.balls[1+i].initSphere(sc*h*Eigen::Vector3d(2.0*std::cos(ang1),2.0*std::sin(ang1), 1), sc*r, moveability); // top ring
    set.conn(1+i, 1 + (i+1)%5) = 2; // around top ring
    set.conn(1+i, 6+i) = 2; // zig
    set.conn(6+i, 1+ (i+1)%5) = 2; // zag

    sc = (6+i)==7 || (6+i) == 8 || (6+i)==9 ? bulge : 1.0;
    if ((6+i) == 8)
      sc *= bulge*bulge;
    moveability = (6+i)==6 || (6+i)==10 ? 0.0 : 1.0;
    set.balls[6+i].initSphere(sc*h*Eigen::Vector3d(2.0*std::cos(ang2),2.0*std::sin(ang2), -1), sc*r, moveability); // bottom ring
    set.conn(6+i, 6 + (i+1)%5) = 2; // around bottom ring
    set.conn(6+i,11) = 2; // bottom fan
  }
  set.balls[11].initSphere(Eigen::Vector3d(0,0,-bulge), bulge*r);

/*  set.balls[0].dest_set = "tree-test";
  set.balls[0].dest_ball_id = 0;
  set.balls[1].dest_set = "hill-test";
  set.balls[1].dest_ball_id = 1;
*/

  set.addLeafBalls({0,1,2,4,5,6,7,8,9,10,11});
  set.leaf_union = false;
  set.render_volume_only = true; 
};


auto hill_testb = [&]()
{
  land.sets.push_back(Landscape::Set("hill-testb", 12));
  Landscape::Set &set = land.sets.back();
  set.colour = Eigen::Vector4d(0.85,0.45,0.25,1);
  const double pi = 3.14159265;

  double h = std::sqrt(1.0 / 5.0);
  double r = std::sqrt(1.0 - h); // 0.743

  // Build hill-test geometry first, then relabel indices so hill-testb is a
  // pure index permutation (TWISTED changes labels only, not shape).
  std::vector<Eigen::Vector3d> base_centres(12, Eigen::Vector3d::Zero());
  std::vector<double> base_radii(12, r);
  std::vector<double> base_mobility(12, 1.0);
  Adj base_conn;
  base_conn.resize(12);

  base_centres[0] = Eigen::Vector3d(0,0,1);
  base_radii[0] = r;
  base_mobility[0] = 0.0;
  for (int i = 0; i < 5; i++)
  {
    float ang1 = (double)i * 2.0*pi/5.0;
    float ang2 = ang1 + pi/5.0;

    base_conn(0, 1+i) = 2; // top fan

    double sc = (1+i)==3 || (1+i)==4 ? bulge : 1.0;
    double moveability = (1+i)==1 || (1+i)==2 || (i+i)==5 ? 0.0 : 1.0;
    if ((1+i)==2)
    {
      double H = std::sqrt(4.0*h*h + (ll-h)*(ll-h));
      double R = std::sqrt(H*H - r*r);
      base_centres[1+i] = (ll/std::sqrt(5.0))*Eigen::Vector3d(2.0*std::cos(ang1),2.0*std::sin(ang1), 1);
      base_radii[1+i] = R;
      base_mobility[1+i] = 0.0;
    }
    else
    {
      base_centres[1+i] = sc*h*Eigen::Vector3d(2.0*std::cos(ang1),2.0*std::sin(ang1), 1);
      base_radii[1+i] = sc*r;
      base_mobility[1+i] = moveability;
    }

    base_conn(1+i, 1 + (i+1)%5) = 2; // around top ring
    base_conn(1+i, 6+i) = 2; // zig
    base_conn(6+i, 1+ (i+1)%5) = 2; // zag

    sc = (6+i)==7 || (6+i) == 8 || (6+i)==9 ? bulge : 1.0;
    if ((6+i) == 8)
      sc *= bulge*bulge;
    moveability = (6+i)==6 || (6+i)==10 ? 0.0 : 1.0;
    base_centres[6+i] = sc*h*Eigen::Vector3d(2.0*std::cos(ang2),2.0*std::sin(ang2), -1);
    base_radii[6+i] = sc*r;
    base_mobility[6+i] = moveability;

    base_conn(6+i, 6 + (i+1)%5) = 2; // around bottom ring
    base_conn(6+i,11) = 2; // bottom fan
  }
  base_centres[11] = Eigen::Vector3d(0,0,-bulge);
  base_radii[11] = bulge*r;
  base_mobility[11] = 1.0;

#if defined(HILL_TESTB_ROT72)
  {
    const double ang = 2.0 * pi / 5.0;
    Eigen::Vector3d axis = Eigen::Vector3d::UnitZ();
#if defined(HILL_TESTB_ROT72_TILTED_AXIS)
    axis = Eigen::Vector3d(1.0, 1.0, 0.5).normalized();
#endif
    Eigen::Matrix3d R = Eigen::AngleAxisd(ang, axis).toRotationMatrix();
    for (auto &c : base_centres)
      c = R * c;
  }
#endif

  // new_index -> old_index in hill-test numbering.
  const int perm[12] = {1, 2, 6, 10, 5, 0, 3, 4, 7, 8, 9, 11};

  for (int ni = 0; ni < 12; ni++)
  {
    int oi = perm[ni];
    set.balls[ni].initSphere(base_centres[oi], base_radii[oi], base_mobility[oi]);
  }
  for (int ni = 0; ni < 12; ni++)
    for (int nj = ni + 1; nj < 12; nj++)
      set.conn(ni, nj) = base_conn(perm[ni], perm[nj]);

  set.addLeafBalls({0,1,2,4,5,6,7,8,9,10,11});
  set.leaf_union = false;
  set.render_volume_only = true; 
};

