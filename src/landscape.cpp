#include "landscape.h"

void Landscape::Set::applyConnectivity()
{
  // adjust spheres to match conn
}

void Landscape::addSetToTypes(Set &set)
{
  for (int i = 0; i<(int)set.balls.size(); i++)
  { 
    Type new_type(set.conn, i);
    bool found = false;
    for (int j = 0; j<types.size(); j++)
    {
      if (types[j].conn == new_type.conn)
      {
        // new_type already exists, so add it in.
        types[j].balls.push_back(&set.balls[i]);
        set.balls[i].type = &types[j];
        found = true;
        break;
      }
    }
    if (!found)
    {
      types.push_back(new_type);
      types.back().balls.push_back(&set.balls[i]);
      set.balls[i].type = &types.back();
    }
  }
}

void Landscape::matchUpDestinationBalls()
{
  for (auto &set: sets)
  {
    for (auto &ball: set.balls)
    {
      if (ball.dest_set == "") // standard recursion
      {
        ball.dest_ball = &ball;
        continue;
      }
      const Type *type = ball.type;
      Set::Ball *dest_ball = nullptr;
      bool found_ball = false;
      for (auto &type_ball: type->balls) // for every ball of this type
      {
        if (type_ball->parent_set->name == ball.dest_set) // check it matches the set name
        {
          Set *par = type_ball->parent_set;
          if (ball.dest_ball_id != -1) // if we have a specific id within the set
          {
            if (&par->balls[ball.dest_ball_id] == type_ball) // then check that id is this ball
            {
              ball.dest_ball = type_ball;
              break;
            }
          }
          else // otherwise just use the first ball that matches the set name
          {
            ball.dest_ball = type_ball;
            break;
          }
        }
      }
    }
  }
}