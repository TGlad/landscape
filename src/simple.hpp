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
    set.balls[2].initSphere(Eigen::Vector3d(-1,0,0), 0.8);
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
    double eq = 0.3; // equatorial distance and radius
    set.balls[0].initSphere(Eigen::Vector3d(0, 0, 1), 1.0);
    set.balls[1].initSphere(Eigen::Vector3d(eq,0, 0), eq);
    set.balls[2].initSphere(Eigen::Vector3d(0,eq, 0), eq);
    set.balls[3].initSphere(Eigen::Vector3d(-eq,0,0), eq);
    set.balls[4].initSphere(Eigen::Vector3d(0,-eq,0), eq);
    set.balls[5].initSphere(Eigen::Vector3d(0, 0,-1), 1.0);
    set.addLeafBall(0,1,2,3);
    set.addLeafBall(2,3,4,5);
    set.leaf_union = true;
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
    set.colour = Eigen::Vector4d(1,0,0,1);
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
    set.balls[0].initSphere(Eigen::Vector3d(0,0,0.35), 0.25);
    set.balls[1].initSphere(Eigen::Vector3d(1,0,0), 1.0);
    set.balls[2].initSphere(Eigen::Vector3d(0,1,0), 1.0);
    set.balls[3].initSphere(Eigen::Vector3d(-1,0,0), 1.0);
    set.balls[4].initSphere(Eigen::Vector3d(0,-1,0), 1.0);
    set.balls[5].initSphere(Eigen::Vector3d(0,0,-0.35), 0.25);

    set.balls[0].dest_set = "tree-tree";
    set.balls[0].dest_ball_id = 0;
    set.addLeafBall(0,1,2,4);
    set.addLeafBall(2,3,4,5);
    set.leaf_union = false;
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
