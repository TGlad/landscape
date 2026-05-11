// ══════════════════════════════════════════════════════════════════
//  Sphere-Inversion Sponge Fractal  —  Shadertoy
//  Algorithm (per DE call):
//    scale = 1
//    repeat 10 times:
//      for each sphere: if p is inside, invert p through that sphere
//                       and accumulate the inversion scale factor
//    return length(p − (1,1,1)) / scale
//
//  Paste the entire file into the Shadertoy "Image" tab and hit ▶.
// ══════════════════════════════════════════════════════════════════

// ─── Fractal parameters ────────────────────────────────────────────
const int   ITERS   = 30;
const bool show_generators = false;//true;
const bool show_arrows = true;

// ─── Ray-march tuning ──────────────────────────────────────────────
const float SKIN        = 0.00125;   // hit (surface) threshold
const int   MAX_STEPS   = 180;     // max ray steps
const float FAR         = 28.0;    // ray kill distance
const float STEP_SCALE  = 0.45;    // DE multiplier < 1 keeps rays safe on fractals
vec4 set_colour = vec4(0,0,0,0);
vec3  ro = vec3(0,0,0);
bool change_colour = false;

float DE_generators(int set, vec3 p)
{
  int child_sets[10] = int[10](set,-1,-1,-1,-1,-1,-1,-1,-1,-1);
  ivec2 child_offsets[10] = ivec2[10](ivec2(0,0),ivec2(0,0),ivec2(0,0),ivec2(0,0),ivec2(0,0),ivec2(0,0),ivec2(0,0),ivec2(0,0),ivec2(0,0),ivec2(0,0));
  int ss[10] = int[10](0,0,0,0,0,0,0,0,0,0);

  int num_siblings[10] = int[10](1,0,0,0,0,0,0,0,0,0);
  int num_child_sets = 1;
  float min_dist = 10000.0;
  int render_set = set;
  for (int s = 0; s<num_child_sets; s++)
  {
    set = child_sets[s];
    ss[set] = s;
    for (int i = SET_OFFSET[set]; i<SET_OFFSET[set] + SET_SIZE[set]; i++)
    {
      int ds = BALLS[i].dest_set;
      if (ds != set)
      {
        // now see if we have to add this
        bool found = false;
        for (int s2 = 0; s2 < num_child_sets; s2++)
        {
          if (child_sets[s2] == ds)
          {
            found = true;
            break;
          }
        }
        if (!found && num_child_sets < 9)
        {
          child_offsets[num_child_sets].x = child_offsets[s].x + 1;
          child_offsets[num_child_sets].y = child_offsets[s].y + num_siblings[child_offsets[num_child_sets].x]++;
          child_sets[num_child_sets++] = ds; 
          ss[ds] = num_child_sets-1;
        }
      }
    }    
  }
  
  for (int s = 0; s<num_child_sets; s++)
  {
    set = child_sets[s];
    for (int i = SET_OFFSET[set]; i<SET_OFFSET[set] + SET_SIZE[set]; i++)
    {
      float len = 4.0;
      float r = 1.0/BALLS[i].curvature;
      vec3 centre = BALLS[i].dir*(BALLS[i].dist + r);
      float offx = len*(float(child_offsets[s].x) - 0.5);
      centre.x -= 0.75*offx;
      centre.y += 0.5*offx;
      centre.z += 2.0*len * (0.5 - (float(child_offsets[s].y)+0.5)/float(num_siblings[child_offsets[s].x])); 
      
      int ds = BALLS[i].dest_set;
      if (ds != set && show_arrows)
      {
        int i2 = BALLS[i].dest_ball;
        int s2 = ss[ds];
        
        float r2 = 1.0/BALLS[i2].curvature;
        vec3 centre2 = BALLS[i2].dir*(BALLS[i2].dist + r2);
        float offx2 = len*(float(child_offsets[s2].x) - 0.5);
        centre2.x -= 0.75*offx2;
        centre2.y += 0.5*offx2;
        centre2.z += 2.0*len * (0.5 - (float(child_offsets[s2].y)+0.5)/float(num_siblings[child_offsets[s2].x])); 
        
        vec3 to2 = centre2 - centre;
        vec3 cent = centre + (ro - p)*0.3; // bring towards camera
        cent += to2*0.05;
        to2 -= to2 * 0.1;
        float d = dot(p - cent, to2) / dot(to2, to2);
        d = max(0.0, min(d, 1.0));
        float dist = length((cent + to2*d) - p) - 0.05;
        vec2 to3 = vec2(0.16, 0.1);
        vec2 ps = vec2(1.0-d, dist);
        float e = max(0.0, min(dot(ps, to3)/dot(to3,to3), 1.0));
        float dist2 = length(to3*e - ps);
        dist = min(dist, dist2);
        if (dist < min_dist)
        {
          min_dist = dist; 
          render_set = -1;
        }
      }
      
      float dist = length(p - centre) - r;
      if (dist < min_dist)
      {
        min_dist = dist;
        render_set = set;
        if (BALLS[i].dest_set != set)
          render_set = BALLS[i].dest_set;
      }
    }
  }
  if (change_colour)
    set_colour = render_set == -1 ? vec4(1,1,1,1) : SET_COLOUR[render_set];
  return min_dist;
}

// ══════════════════════════════════════════════════════════════════
//  Distance Estimator
//  For each iteration we walk through all 16 inversion spheres.
//  Whenever p lies inside a sphere we reflect/invert it outward and
//  multiply `scale` by the same dilation factor.
// ══════════════════════════════════════════════════════════════════
float DE(vec3 p) 
{
  int set = 0;
  if (show_generators)
    return DE_generators(set, p);
  float scale = 1.0;

  float rad0 = SKIN;
  float rad = SKIN;
  int location_ids[MAX_BALLS_PER_SET];
  int loc_balls[MAX_BALLS_PER_SET];
  int num_locations = 0;
  for (int i = 0; i<SET_SIZE[set]; i++)
  {
    int j = SET_OFFSET[set] + i;
    location_ids[i] = LOCATION_OFFSETS[j];
    if (location_ids[i] < LOCATION_OFFSETS[j + 1])
      loc_balls[num_locations++] = i;
  }

  int dest_set = -1; // set;
  int dest_offset = 0;
  int dest_num = 0;
  int dest_i = 0;
  for (int n = 0; n < ITERS; n++) 
  {
    bool found = false;

    // substitution sets
    // TODO: perhaps generate substitution set offset and set_offset separately to avoid continue below
    for (int i = SET_OFFSET[set]; i<SET_OFFSET[set] + SET_SIZE[set]; i++)
    {
      if (BALL_MOBIUS[i] == -1)
        continue;
      int I = i - SET_OFFSET[set];

      float r = 1.0/BALLS[i].curvature;
      vec3 centre = BALLS[i].dir*(BALLS[i].dist + r);
      vec3 off = p - centre;
      float d2 = dot(off,off);
      if (d2 < 4.0*rad*rad) // TODO: the 4x works better and is faster, but why?
      {
        found = false;
        break;
      }
      if (d2 < r*r)
      {
        // if it has a dest then invert as normal, 
        // specify the dest_set, but don't change set until you exit A,B
        // if you exit A,B then swap set
        // if you get to tree (OVERLAPS[0].ball_0) then swap to overlap set instead 
        if (location_ids[I] >= LOCATION_OFFSETS[i+1])
        {
          found = false;
          for (int j = OVERLAP_OFFSET[set]; j<OVERLAP_OFFSET[set+1]; j++)
          {
            int b0 = OVERLAPS[j].ball_0;
            int ds_0 = BALLS[b0].dest_set;
            int b1 = OVERLAPS[j].ball_1;
            int ds_1 = BALLS[b1].dest_set;
            if (dest_set != ds_0 && dest_set != ds_1)
              continue;
              
            int bi = dest_set == ds_0 ? b1 : b0;
      
            float r0 = 1.0/BALLS[bi].curvature;
            vec3 centre0 = BALLS[bi].dir*(BALLS[bi].dist + r0);
            vec3 off0 = p - centre0;
            float d0 = dot(off0,off0);

            if (d0 < r0*r0) // inside both spheres
            {
              set = OVERLAPS[j].dest_set;
              dest_set = set;
              num_locations = 0;
              // don't do Mobius transform quite yet
              for (int i = 0; i<SET_SIZE[set]; i++)
              {
                int j = SET_OFFSET[set] + i;
                location_ids[i] = LOCATION_OFFSETS[j];
                if (location_ids[i] < LOCATION_OFFSETS[j + 1])
                  loc_balls[num_locations++] = i;
              }
              found = true;
              break; // just to be sure
            }
          }
          if (found)
            break;
        }
        
        // check if we are switching
        found = true;

        p = off;
        float len = sqrt(d2);
        float lmin = r*r / (len + rad);
        float lmax = r*r / (len - rad);
        rad = (lmax - lmin)/2.0;
        p = p * (lmax + lmin)/(2.0*len) + centre;


        if (location_ids[I] >= LOCATION_OFFSETS[i+1])
        {
          dest_set = BALLS[i].dest_set;
          dest_i = i;
          dest_offset = BALLS[i].neighbours_offset;
          dest_num = BALLS[i].num_neighbours;        
          break;
        }

        for (int k = 0; k<num_locations; k++)
        {
          int j = loc_balls[k];
          if (location_ids[j] < NUM_LOCATIONS && LOCATIONS[location_ids[j]] == I)
            location_ids[j]++;
        }
        break;
      }
    }
    if (found)
      continue;
      
    // standard spheres
    for (int i = SET_OFFSET[set]; i<SET_OFFSET[set] + SET_SIZE[set]; i++)
    {
      if (BALL_MOBIUS[i] != -1)
        continue;
      int I = i - SET_OFFSET[set];

      float r = 1.0/BALLS[i].curvature;
      vec3 centre = BALLS[i].dir*(BALLS[i].dist + r);
      vec3 off = p - centre;
      float d2 = dot(off,off);
      if (d2 < 4.0*rad*rad) // TODO: the 4x works better and is faster, but why?
      {
        found = false;
        break;
      }
      if (d2 < r*r)
      {
        // check if we are switching
        found = true;
        if (set != dest_set)
        {
          bool inside = i == dest_i;
          for (int k = 0; k<BALLS[dest_i].num_neighbours; k++)
          {
            int n_id = NEIGHBOURS[BALLS[dest_i].neighbours_offset + k];
            if (i == n_id)
              inside = true;
          }
          if (dest_set != -1 && !inside)
          {
            p = applyMobius(BALL_MOBIUS[dest_i], p, rad);
            set = dest_set;
            num_locations = 0;
            for (int i = 0; i<SET_SIZE[set]; i++)
            {
              int j = SET_OFFSET[set] + i;
              location_ids[i] = LOCATION_OFFSETS[j];
              if (location_ids[i] < LOCATION_OFFSETS[j + 1])
                loc_balls[num_locations++] = i;
            }
            break;              
          }
        }
        p = off;
        float len = sqrt(d2);
        float lmin = r*r / (len + rad);
        float lmax = r*r / (len - rad);
        rad = (lmax - lmin)/2.0;
        p = p * (lmax + lmin)/(2.0*len) + centre;

        for (int k = 0; k<num_locations; k++)
        {
          int j = loc_balls[k];
          if (location_ids[j] < NUM_LOCATIONS && LOCATIONS[location_ids[j]] == I)
            location_ids[j]++;
        }
        break;
      }
    }
    if (!found)
      break;
  }

  if (change_colour)
    set_colour = SET_COLOUR[set];

  float dist = 1000.0;
  if (LEAF_SIZE[set] > 0)
  {
    dist = LEAF_UNION[set] ? 1000.0 : -1000.0;
    for (int i = LEAF_OFFSET[set]; i<LEAF_OFFSET[set]+LEAF_SIZE[set]; i++)
    {
      float r = 1.0/LEAF_BALLS[i].curvature;
      vec3 centre = LEAF_BALLS[i].dir*(LEAF_BALLS[i].dist + r);
      if (LEAF_UNION[set])
        dist = min(dist, length(p - centre)-r);
      else // intersection
        dist = max(dist, length(p - centre)-r);
    }
    dist += rad; // to get rid of the skin width on volumes
  }
/*    float num = 0.0;
    for (int i = SET_OFFSET[set]; i<SET_OFFSET[set] + SET_SIZE[set]; i++)
    {
      if (BALLS[i].curvature != 0.0)
      {
        float r = 1.0/BALLS[i].curvature;
        vec3 centre = BALLS[i].dir*(BALLS[i].dist + r);
        dist = min(dist, length(p - centre) - 0.1);
      }
    }*/
  if (!VOLUME_ONLY[set])
  {
    vec3 mid = vec3(0,0,0);
    float num = 0.0;
    for (int i = SET_OFFSET[set]; i<SET_OFFSET[set] + SET_SIZE[set]; i++)
    {
      if (BALLS[i].curvature != 0.0)
      {
        float r = 1.0/BALLS[i].curvature;
        mid += BALLS[i].dir*(BALLS[i].dist + r);
        num++;
      }
    }
    mid /= num;
    dist = min(dist, length(p - mid));
  }
  return dist * rad0/rad;
}


// ══════════════════════════════════════════════════════════════════
//  Central-difference surface normal
// ══════════════════════════════════════════════════════════════════
vec3 calcNormal(vec3 p) {
    const vec2 e = vec2(SKIN*0.5, 0.0);
    return normalize(vec3(
        DE(p + e.xyy) - DE(p - e.xyy),
        DE(p + e.yxy) - DE(p - e.yxy),
        DE(p + e.yyx) - DE(p - e.yyx)
    ));
}


// ══════════════════════════════════════════════════════════════════
//  Approximate ambient occlusion
//  March short distances along the normal; deficit ↔ occlusion.
// ══════════════════════════════════════════════════════════════════
float calcAO(vec3 p, vec3 n) {
    float occ = 0.0, w = 1.0;
    for (int i = 1; i <= 5; i++) {
        float h = 0.03 * float(i);
        occ += w * max(h - DE(p + n * h), 0.0);
        w   *= 0.5;
    }
    return clamp(1.0 - 5.0 * occ, 0.0, 1.0);
}


// ══════════════════════════════════════════════════════════════════
//  Ray marcher
//  Returns hit distance t, or −1 on miss.
// ══════════════════════════════════════════════════════════════════
float rayMarch(vec3 ro, vec3 rd) {
    float t = 0.05;
    change_colour = true;
    for (int i = 0; i < MAX_STEPS; i++) {
        float d = DE(ro + rd * t);
        if (d < SKIN)  return t + d-SKIN;
        if (t > FAR)   return -1.0;
        t += d * STEP_SCALE;
    }
    return -1.0;
}


// ══════════════════════════════════════════════════════════════════
//  Shadertoy entry point
// ══════════════════════════════════════════════════════════════════
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // ── Screen-space UV, Y-up, aspect-corrected ───────────────────
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;

    // ── Slowly orbiting camera ─────────────────────────────────────
    float a  = iTime * 0.20;
    ro = vec3(sin(a) * 5.5,
                    cos(a) * 5.5,
                    show_generators ? 0.0 : 1.8 + sin(a * 0.37) * 1.4);
    vec3  ta = vec3(0.0,0.0,0.0);                         // look-at target

    vec3  fwd = normalize(ta - ro);
    vec3  rgt = normalize(cross(fwd, vec3(0.0, 0.0, 1.0)));
    vec3  upv = cross(rgt, fwd);
    float fl  = show_generators ? 0.6 : 3.8;                            // focal length / zoom
    vec3  rd  = normalize(fl * fwd + uv.x * rgt + uv.y * upv);

    // ── Sky gradient ──────────────────────────────────────────────
    vec3 sky = mix(
        vec3(0.04, 0.06, 0.14),
        vec3(0.25, 0.40, 0.62),
        clamp(0.5 + 0.7 * rd.z, 0.0, 1.0)
    );
    vec3 col = sky;

    // ── Ray march ─────────────────────────────────────────────────
    float t = rayMarch(ro, rd);
    change_colour = false;

    if (t > 0.0) {
        vec3  p  = ro + rd * t;
        vec3  n  = calcNormal(p);
        float oc = calcAO(p, n);

        // Diffuse + specular + Fresnel rim
        vec3  ldir = normalize(vec3(0.7, 0.5, 1.2));
        float diff = max(dot(n, ldir), 0.0);
        float spec = pow(max(dot(reflect(-ldir, n), -rd), 0.0), 48.0);
        float fres = pow(1.0 - max(dot(-rd, n), 0.0), 4.0);

        col  = set_colour.xyz * (0.08 * oc + 0.92 * diff * oc);
        col += vec3(0.95, 0.90, 0.80) * spec * 0.55;
        col  = mix(col, vec3(0.75, 0.88, 1.00), fres * 0.18);

        // Atmospheric depth fog
      //  col = mix(col, sky, 1.0 - exp(-0.032 * t * t));
    }

    // ── Gamma correction (sRGB) ────────────────────────────────────
    col = pow(clamp(col, 0.0, 1.0), vec3(0.4545));

    fragColor = vec4(col, 1.0);
}
