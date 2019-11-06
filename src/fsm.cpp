#include <iostream>
#include <cassert>
#include "fsm.h"
#include "helpers.h"

using nlohmann::json;
using std::string;
using std::vector;


/**
 * Helper function to print Movement structure data
 * @param m The movement to print
 */
static void printMovement(const Movement & m)
{
  double tot_v = sqrt(pow(m.vx, 2) + pow(m.vy, 2));
  double tot_a = sqrt(pow(m.ax, 2) + pow(m.ay, 2));
  std::cout << "Movement. yaw: " << m.yaw << " x: " << m.x << " y: " << m.y << " vx: " << m.vx <<
      " vy: " << m.vy << " tot_v: " << tot_v << " ax: " << m.ax << " ay: " << m.ay << " tot_a: " <<
      tot_a << std::endl;
}

/**
 * Helper function to translate World coordinates to local car coordinates
 * @param x The world X coordinate
 * @param y The world Y coordinate
 * @param ref_yaw The reference point world YAW for translation between car and world coords.
 * @param ref_x The reference point world X for translation between car and world coords.
 * @param ref_y The reference point world Y for translation between car and world coords.
 * @return A vector of X, Y local car coords.
 */
inline static std::vector<double>
translateToCarCoords(double x, double y, double ref_yaw, double ref_x, double ref_y)
{
  double shift_x = x - ref_x;
  double shift_y = y - ref_y;
  double nx = (shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw));
  double ny = (shift_x * sin(-ref_yaw) + shift_y * cos(-ref_yaw));
  return {nx, ny};
}

/**
 * Helper function to translate car local coordinates to world coords
 * @param x The car local X coordinate
 * @param y The car local Y coordinate
 * @param ref_yaw The reference point world YAW for translation between car and world coords.
 * @param ref_x The reference point world X for translation between car and world coords.
 * @param ref_y The reference point world Y for translation between car and world coords.
 * @return A vector of X, Y world coords.
 */
static std::vector<double> translateToWorldCoords(double x, double y, double ref_yaw, double ref_x, double ref_y)
{
  double x_ref = x;
  double y_ref = y;
  // Re-translate to original world coordinates
  x = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
  y = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

  x += ref_x;
  y += ref_y;
  return {x, y};
}

/**
 * Constructor for the FSM
 * @param map The Map structure containing the waypoint world maps of the track
 */
FSM::FSM(const Map & map) : 
  state_(nullptr), cxt_(nullptr), cur_policy_(KEEP_LANE), map_(map),
  cur_trg_speed_(0), first_time_(true)
{
  // Make sure compilation fails if we try to use an undefined starting lane
  static_assert(STARTING_LANE == LEFT_LANE_D or STARTING_LANE == CENTER_LANE_D or STARTING_LANE == RIGHT_LANE_D);
  // Set current lane based on STARTING_LANE constexpr
  if (STARTING_LANE == LEFT_LANE_D) lane_ = LEFT;
  else if (STARTING_LANE == CENTER_LANE_D) lane_ = CENTER;
  else lane_ = RIGHT;
}

/**
 * Method to compute proper steps to match a given speed, without violating any constraints
 * @param trg_v The target speed
 */
void MatchSpeed::matchTargetSpeed(double trg_v)
{
  if (DEBUG) std::cout << "Target speed to match: " << trg_v << std::endl;
  forward_traj_.clear();
  cur_traj_idx_ = 0;

  Movement prev_m = fsm_.getLastPlannedMovement();

  double x_step = 0;
  double vx = prev_m.vx;
  double vy = prev_m.vy;
  double ax = prev_m.ax;
  double ay = prev_m.ay;
  double tot_v = sqrt(pow(vx, 2) + pow(vy, 2));
  double tot_a = sqrt(pow(ax, 2) + pow(ay, 2));

  double delta_v = trg_v - tot_v;
  // If we are already at the target speed, with zero acceleration, then nothing to do
  if (delta_v > 0 and delta_v < 0.1 and isZero(tot_a)) return;

  const bool accel = delta_v >= 0;

  if (DEBUG) std::cout << "MatchTargetSpeed: initial speed: " << tot_v << " initial acc: " << tot_a << std::endl;

  // The flip point: speed at which to switch from accelerating to decelerating or vice versa
  double flip_point = tot_v + delta_v / 2;
  if (DEBUG) std::cout << "MatchTargetSpeed: initial flip_point " << flip_point << std::endl;

  // How much do we need to shed to go from max_acc to the initial_one?
  double delta_speed_shed_to_init = 0;
  {
    int shed_acc_steps = ceil(((MAX_ACC - tot_a) / MAX_JERK) / TIME_RES);
    double a_shed = 0;
    for (int i = 0; i < shed_acc_steps; ++i) {
      delta_speed_shed_to_init += a_shed * TIME_RES;
      a_shed += MAX_JERK * TIME_RES;
    }
  }

  double cur_x = prev_m.x;
  double cur_y = prev_m.y;
  double cur_speed = tot_v;
  double prev_speed = cur_speed;
  double cur_acc = tot_a;
  double cur_yaw = prev_m.yaw;
  double prev_x = cur_x;
  double prev_y = cur_y;
  if (DEBUG) std::cout << "MatchTargetSpeed: " << (accel? "accelerating": "decelerating") << std::endl;

  while (accel? trg_v - cur_speed > 0: cur_speed - trg_v > 0) {
    cur_x       += vx * TIME_RES;
    cur_y       += vy * TIME_RES;
    vx          += ax * TIME_RES;
    vy          += ay * TIME_RES;
    cur_speed   = sqrt(pow(vx, 2) + pow(vy, 2));

    bool flipped = (accel)? cur_speed > flip_point: cur_speed < flip_point;

    // Break condition, if we came close to our target
    if (flipped and (accel? cur_speed < prev_speed: prev_speed < cur_speed)) {
      if (DEBUG) std::cout << "MatchTargetSpeed: cur_speed inversion at " << cur_speed << "\n";
      break;
    }

    double jerk  = (accel)? MAX_JERK: -MAX_JERK;
    jerk         = (flipped)? -jerk: jerk;
    double jx    = jerk * cos(cur_yaw);
    double jy    = jerk * sin(cur_yaw);
    double n_ax  = ax + jx * TIME_RES;
    double n_ay  = ay + jy * TIME_RES;
    cur_acc      = sqrt(pow(n_ax, 2) + pow(n_ay, 2));
    if (cur_acc >= MAX_ACC) {
      cur_acc = MAX_ACC;
      ax = cur_acc * cos(cur_yaw);
      ay = cur_acc * sin(cur_yaw);
      // Update flip point as well
      flip_point = trg_v + (accel? -1: 1) * delta_speed_shed_to_init;
      if (DEBUG) std::cout << "MatchTargetSpeed: updated flip_point " << flip_point << std::endl;
    } else {
      ax = n_ax;
      ay = n_ay;
    }
    cur_yaw = atan2(cur_y - prev_y, cur_x - prev_x);
    prev_x = cur_x;
    prev_y = cur_y;

    forward_traj_.emplace_back();
    TrajNode & tn = forward_traj_.back();
    tn.x = cur_x;
    tn.y = cur_y;

    if (DEBUG) {
      std::cout << "Cur speed: " << cur_speed << " Cur acc: " << cur_acc << " Cur x: " << cur_x <<
          " y: " << cur_y << " vx: " << vx << " vy: " << vy << " ax: " << ax << " ay: " << ay <<
          " flipped: " << flipped << std::endl;
    }
    prev_speed = cur_speed;
  }
  if (DEBUG) std::cout << "MatchTargetSpeed: computed trajectory of " << forward_traj_.size() << " steps\n";
}



/**
 * Function to return the closest car ahead or behind, in the lane defined by d_offset
 * @param p The most recently planned position data.
 * @param cxt The context containing the sensor fusion data
 * @param ahead Whether we should be looking for cars ahead of us (True) or behind us (False)
 * @param d_offset The offset in D we should be looking at, relative to current position (typically -4, 0 or +4)
 * @return The sensor data containing information concerning the closest car ahead/behind in the request d_offset lane
 */
static SensorData getSensorFusion(const PositionData & p, Context & cxt, bool ahead, double d_offset)
{
  // If we are looking at the side lanes, look a bit further ahead (not behind!)
  double factor = 1.0;
  if (d_offset > HALF_LANE_WIDTH or d_offset < -HALF_LANE_WIDTH) factor = 1.5;

  // Get data for car closest to us in s, going forward
  int min_idx  = -1;
  double min_diff_s = MAX_S;
  for (int i = 0; i < cxt.sensor_fusion.size(); ++i) {
    // Only interested about cars in the same lane for the given offset
    double car_d = cxt.sensor_fusion[i][6]; // get the d value
    if (car_d < p.d + d_offset - HALF_LANE_WIDTH or
        car_d > p.d + d_offset + HALF_LANE_WIDTH) continue;

    double own_s  = cxt.s;
    double car_s  = cxt.sensor_fusion[i][5]; // get the s value

    // Handle wraparound
    if (ahead) {
      if (car_s  < own_s) car_s += MAX_S;
    } else {
      if (car_s  > own_s) own_s += MAX_S;
    }
    
    // Only interested about cars in a reasonable distance
    double diff_s = fabs(car_s - own_s);
    double margin = SENSOR_MAX_DIST;
    if (car_s - own_s > 0) margin *= factor;
    if (diff_s > margin) continue;

    if (diff_s < min_diff_s) {
      min_diff_s  = diff_s;
      min_idx     = i;
    }
  }
  SensorData sd; 
  sd.diff_s = min_diff_s;
  // Are there any cars detectable ahead/behind us?
  if (min_idx == -1) {
    sd.s = -1;
    sd.d = p.d;
    return sd;
  }
  sd.id = cxt.sensor_fusion[min_idx][0];
  sd.x  = cxt.sensor_fusion[min_idx][1];
  sd.y  = cxt.sensor_fusion[min_idx][2];
  sd.vx = cxt.sensor_fusion[min_idx][3];
  sd.vy = cxt.sensor_fusion[min_idx][4];
  sd.v  = sqrt(pow(sd.vx, 2) + pow(sd.vy, 2));
  sd.s  = cxt.sensor_fusion[min_idx][5];
  sd.d  = cxt.sensor_fusion[min_idx][6];
  return sd;
}

/**
 * Function that decides whether we can change lane, based on cars ahead and behind
 * @param ahead The sensor data of the closest car ahead of us
 * @param behind The sensor data of the closest car behind us
 * @return True if we can change lane, False otherwise
 */
bool FSM::canChangeLane(SensorData & ahead, SensorData & behind)
{
  double time_before_change = cxt_->path_x.size() * TIME_RES;
  double time_to_change = LANE_CHANGE_FCT / cur_trg_speed_; // time to change lane with current speed
  // First check the car behind
  {
    if (behind.diff_s < LANE_CHANGE_MAR) return false;
    double dv = behind.v - cur_trg_speed_;
    double future_diff_s = behind.diff_s - LANE_CHANGE_MAR - time_before_change * dv;
    if (future_diff_s < dv * time_to_change) return false;
  }
  // Now check the car ahead
  {
    if (ahead.diff_s < LANE_CHANGE_MAR) return false;
    double dv = cur_trg_speed_ - ahead.v;
    double future_diff_s = ahead.diff_s - LANE_CHANGE_MAR - time_before_change * dv;
    if (future_diff_s < dv * time_to_change) return false;
  }
  return true;
}

/**
 * Helper function for printing out the current lane
 * @param lane The lane to print
 * @return lane string to print
 */
static std::string printLane(Lane lane)
{
  switch (lane) {
  case CENTER:
    return "CENTER";
  case LEFT:
    return "LEFT";
  case RIGHT:
    return "RIGHT";
  default:
    assert(false);
  }
  return "";
}

/**
 * Helper function for printing the current Policy
 * @param policy The policy to print
 * @return policy string to print
 */
static std::string printPolicy(Policy policy)
{
  switch (policy) {
  case KEEP_LANE:
    return "KEEP_LANE";
  case CHANGE_LANE:
    return "CHANGE_LANE";
  case MATCH_SPEED:
    return "MATCH_SPEED";
  default:
    assert(false);
  }
  return "";
}

/**
 * Function to check which of the available lanes is best to use
 * @param p The most recently planned position data (end of path_x/y or next_x/y_vals) vectors
 * @return The sensor data to use to determine which state to switch to. (if sd.s == -1 then there is no car present)
 */
SensorData FSM::checkLanes(const PositionData & p)
{
  if (DEBUG) std::cout << "checkLanes: current lane: " << printLane(lane_) << std::endl;
  // Check same lane
  SensorData best = getSensorFusion(p, *cxt_, true, 0); best.d = p.d; best.policy = KEEP_LANE;

  if (best.s == -1)  {
    best.v = MAX_SPEED * MAX_SPEED_MAR;
    // If lane is not center, try to see if we can move to the center one
    if (lane_ == CENTER) {
      if (DEBUG) std::cout << "Staying in the same lane, since there is no car there\n";
      return best;
    }
  }

  if (DEBUG and best.s != -1)
      std::cout << "Car in the same lane: " << best.id << " s: " << best.s << " d: " << best.d << " v: " << best.v << std::endl;

  // Check lane to our right
  if (lane_ == LEFT or lane_ == CENTER) {
    auto right_ahead = getSensorFusion(p, *cxt_, true,  LANE_WIDTH);
    auto right_backw = getSensorFusion(p, *cxt_, false, LANE_WIDTH);
    // Check if there is enough space
    if (canChangeLane(right_ahead, right_backw)) {
      // If there are no cars present ahead right or if their speed is higher than the best, chose right
      if (right_ahead.s == -1 or best.v < 0.95 * right_ahead.v) {
        best = right_ahead;
        best.policy = CHANGE_LANE;
        best.d = (lane_ == LEFT)? CENTER_LANE_D: RIGHT_LANE_D;
        if (best.s == -1)  {
          if (DEBUG) std::cout << "Moving to right lane, since there is no car there\n";
          best.v = cur_trg_speed_;
          return best;
        }
        if (DEBUG) std::cout << "Car in our right lane: " << best.id << " s: " << best.s << " d: " << best.d << 
                     " v: " << best.v << std::endl;
      }
    }
  }

  // Check lane to our left
  if (lane_ == RIGHT or lane_ == CENTER) {
    auto left_ahead = getSensorFusion(p, *cxt_, true,  -LANE_WIDTH);
    auto left_backw = getSensorFusion(p, *cxt_, false, -LANE_WIDTH);
    // Check if there is enough space
    if (canChangeLane(left_ahead, left_backw)) {
      // If there are no cars present ahead left or if their speed is higher than the best, chose left
      if (left_ahead.s == -1 or best.v < 0.95 * left_ahead.v) {
        best = left_ahead;
        best.policy = CHANGE_LANE;
        best.d = (lane_ == RIGHT)? CENTER_LANE_D: LEFT_LANE_D;
        if (best.s == -1)  {
          if (DEBUG) std::cout << "Moving to left lane, since there is no car there\n";
          best.v = cur_trg_speed_;
          return best;
        }
        if (DEBUG) std::cout << "Car in our left lane: " << best.id << " s: " << best.s << " d: " << best.d << 
                     " v: " << best.v << std::endl;
      }
    }
  }

  if (DEBUG) {
    std::cout << "Applying policy: " << printPolicy(best.policy) << std::endl;
    std::cout << "  Target lane with d: " << best.d << std::endl;
  }

  return best;
}

/**
 * Method to update the internal lane state
 * @param d The new D value of the current position
 */
void FSM::updateLane(double d)
{
  assert(d >= 0 and d <= 3 * LANE_WIDTH);
  if (d >= 0 and d < LANE_WIDTH) lane_ = LEFT;
  else if (d >= LANE_WIDTH and d < 2 * LANE_WIDTH) lane_ = CENTER;
  else if (d >= 2 * LANE_WIDTH) lane_ = RIGHT;
}

/**
 * Update method to be called for every update of the FSM
 */
void State::update()
{
  auto p = fsm_.getMostRecentNotProcessed();
  // Clear old processed data on the FSM
  fsm_.clearProcessed();
  // Apply the new trajectory by creating the new next_x_vals and next_y_vals
  fsm_.applyTrajectory(p, forward_traj_, cur_traj_idx_);
}

/**
 * Change Lane STATE
 * @param fsm The FSM object
 * @param trg_speed The target speed to use
 * @param sd The sensor data for any car we detected, or at least one showing the lane to change to (should not be null)
 */
ChangeLane::ChangeLane(FSM & fsm, double trg_speed, SensorData * sd) : State(fsm, trg_speed)
{
  if (DEBUG) std::cout << "ChangeLane: creating new State\n";
  fsm_.cur_policy_ = CHANGE_LANE;
  auto p = fsm_.getMostRecentNotProcessed();
  fsm_.clearProcessed();
  fsm_.cur_trg_speed_ = trg_speed_;
  p.d = sd->d;
  fsm_.applyTrajectory(p, forward_traj_, cur_traj_idx_);
}

/**
 * Method determining next policy
 * @return the next policy to follow
 */
Policy ChangeLane::nextPolicy()
{
  if (cur_traj_idx_ >= forward_traj_.size() + 2 * TRAJ_STEPS) {
    if (DEBUG) std::cout << "ChangeLane: completed maneuver. Next policy is KEEP_LANE\n";
    return KEEP_LANE;
  }
  if (DEBUG) std::cout << "ChangeLane: next policy remains CHANGE_LANE\n";
  return CHANGE_LANE;
}

/**
 * Match Speed STATE
 * @param fsm The FSM object
 * @param trg_speed The target speed to match
 * @param sd Optionally the sensor data for any car we detected (nullptr otherwise)
 */
MatchSpeed::MatchSpeed(FSM & fsm, double trg_speed, SensorData * sd) : State(fsm, trg_speed)
{
  if (DEBUG) std::cout << "MatchSpeed: creating new State\n";
  fsm_.cur_policy_ = MATCH_SPEED;
  auto p = fsm_.getMostRecentNotProcessed();
  fsm_.clearProcessed();
  matchTargetSpeed(trg_speed_);
  fsm_.cur_trg_speed_ = trg_speed_;
  // Apply the new trajectory by creating the new next_x_vals and next_y_vals
  fsm_.applyTrajectory(p, forward_traj_, cur_traj_idx_);
}

/**
 * Method determining next policy
 * @return the next policy to follow
 */
Policy MatchSpeed::nextPolicy()
{
  if (cur_traj_idx_ >= forward_traj_.size() + 2 * TRAJ_STEPS) {
    if (DEBUG) std::cout << "MatchSpeed: completed maneuver. Next policy is KEEP_LANE\n";
    return KEEP_LANE;
  }
  if (DEBUG) std::cout << "MatchSpeed: next policy remains MATCH_SPEED\n";
  return MATCH_SPEED;
}


//
//
// FSM implementation
//
//

/**
 * Method to get the last movement that was planned
 * @return Movement structure containing last movement
 */
Movement FSM::getLastPlannedMovement()
{
  size_t sz = next_x_vals_.size();
  Movement m;
  if (sz < 3) {
    m.yaw = deg2rad(cxt_->yaw);
    m.x = sz? next_x_vals_[sz-1]: cxt_->x;
    m.y = sz? next_y_vals_[sz-1]: cxt_->y;
    m.vx = m.vy = 0;
    m.ax = m.ay = 0;
  } else {
    double x1 = next_x_vals_[sz-1];
    double y1 = next_y_vals_[sz-1];
    double x2 = next_x_vals_[sz-2];
    double y2 = next_y_vals_[sz-2];
    double x3 = next_x_vals_[sz-3];
    double y3 = next_y_vals_[sz-3];
    // TODO: assuming that is x1 == x2 then we are at a right angle to the x axis.
    if (x1 == x2) {
      m.yaw = cxt_->yaw;
      if (DEBUG) std::cout << "ATAN2 for same X's not defined, using: " << m.yaw << std::endl;
    } else {
      m.yaw = atan2(y1 - y2, x1 - x2);
    }
    m.x = x1;
    m.y = y1;
    m.vx = (x1 - x2) / TIME_RES;
    m.vy = (y1 - y2) / TIME_RES;
    double prev_vx = (x2 - x3) / TIME_RES;
    double prev_vy = (y2 - y3) / TIME_RES;
    m.ax = (m.vx - prev_vx) / TIME_RES;
    m.ay = (m.vy - prev_vy) / TIME_RES;
  }
  return m;
}

/**
 * Method to generate a spline based on current position, previous couple ones,
 * and certain points on the road ahead, spaced by 30 meters
 * @param d. The target d position for the car
 * @param prev_m. The previous movement details
 * @param ref_x. Will output the x for translation between world/car coords
 * @param ref_y. Will output the y for translation between world/car coords
 * @param ref_yaw. Will output the ref_yaw for translation between world/car coords
 * @return Will return a spline created for use with local car coords
 */
tk::spline * FSM::createLocalSpline(double d, const Movement & prev_m,
    double & ref_x, double & ref_y, double & ref_yaw)
{
  Context & cxt = *cxt_;
  size_t sz = next_x_vals_.size();

  // Here we need to start with two simple points close to our reference location
  vector<double> xs, ys;
  ref_yaw = deg2rad(cxt.yaw);
  ref_x   = cxt.x;
  ref_y   = cxt.y;

  if (sz < 2) {
    double prev_x = ref_x - cos(ref_yaw);
    double prev_y = ref_y - sin(ref_yaw);
    xs.push_back(prev_x); ys.push_back(prev_y);
  } else {
    ref_x = next_x_vals_[sz-1];
    ref_y = next_y_vals_[sz-1];
    double ref_x_prev = next_x_vals_[sz-2];
    double ref_y_prev = next_y_vals_[sz-2];
    ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);
    xs.push_back(ref_x_prev); ys.push_back(ref_y_prev);
  }
  xs.push_back(ref_x); ys.push_back(ref_y);

  auto cur_sd = getFrenet(ref_x, ref_y, ref_yaw, map_.x, map_.y);
  double cur_d = cur_sd[1];

  double d1, d2, d3;
  d1 = d2 = d3 = d;
  double s1, s2 , s3;
  // Just picking some waypoints further ahead
  s1 = cur_sd[0] + 30;
  s2 = s1 + 30;
  s3 = s2 + 30;
  // Now we will add some more points further away
  auto next_wp1 = getXY(s1, d1, map_.s, map_.x, map_.y);
  auto next_wp2 = getXY(s2, d2, map_.s, map_.x, map_.y);
  auto next_wp3 = getXY(s3, d3, map_.s, map_.x, map_.y);
  xs.push_back(next_wp1[0]); ys.push_back(next_wp1[1]);
  xs.push_back(next_wp2[0]); ys.push_back(next_wp2[1]);
  xs.push_back(next_wp3[0]); ys.push_back(next_wp3[1]);

  // Now we need to transform these points to the car-local coordinates:
  // The reference point is at (0,0) and the x axis points to our yaw angle.
  // This way we can use the spline as f(x), and it would still behave as a function (single value of y for x)
  for (int i = 0; i < xs.size(); ++i) {
    auto car_xy = translateToCarCoords(xs[i], ys[i], ref_yaw, ref_x, ref_y);
    xs[i] = car_xy[0];
    ys[i] = car_xy[1];
  }

  tk::spline * spl = new tk::spline();
  spl->set_points(xs, ys);
  return spl;
}

/**
 * Method to populate the future X,Y trajectory based on a Spline interpolation
 * @param d Current d Frenet value
 * @param ft The forward trajectory, if one is computed (empty otherwise)
 * @param cti The current trajectory index. Will be incremented by the steps taken.
 */
void FSM::createXYFromSpline(double d, vector<TrajNode> & ft, size_t & cti)
{
  Movement prev_m = getLastPlannedMovement();
  auto init_sd = getFrenet(prev_m.x, prev_m.y, prev_m.yaw, map_.x, map_.y);
  if (DEBUG) std::cout << "Target speed for Spline: " << cur_trg_speed_ << std::endl;

  size_t sz = next_x_vals_.size();

  // Create a local spline for use with local car coords
  double ref_x, ref_y, ref_yaw;
  std::unique_ptr<tk::spline> spl(createLocalSpline(d, prev_m, ref_x, ref_y, ref_yaw));

  // How many point do I need to sample from the spline so that we still keep under our
  // velocity and acceleration constraints?
  double trg_x = LANE_CHANGE_FCT;
  double trg_y = (*spl)(trg_x);
  double trg_dist = sqrt(pow(trg_x, 2) + pow(trg_y, 2));
  double N = trg_dist / (TIME_RES * cur_trg_speed_);
  double x_add_on = 0; // Always in car-local coords
  double prev_x_point = prev_m.x;
  double prev_y_point = prev_m.y;
  double cur_yaw = ref_yaw;
  for (size_t t = sz; t < TRAJ_STEPS; ++t) {
    double x_point = x_add_on +  trg_x / N;
    // Translate forward trajectory to car local coords, if present
    if (cti < ft.size()) {
      const TrajNode & tn = ft[cti];
      auto car_xy = translateToCarCoords(tn.x, tn.y, ref_yaw, ref_x, ref_y);
      x_point = car_xy[0];
    };
    ++cti;

    double y_point = (*spl)(x_point);
    auto world_xy = translateToWorldCoords(x_point, y_point, ref_yaw, ref_x, ref_y);
    x_add_on = x_point;
    x_point = world_xy[0];
    y_point = world_xy[1];

    // Update cur_yaw
    cur_yaw = atan2(y_point - prev_y_point, x_point - prev_x_point);
    prev_x_point = x_point;
    prev_y_point = y_point;

    // Populate the rest of the values
    next_x_vals_.push_back(x_point);
    next_y_vals_.push_back(y_point);
    auto cur_sd = getFrenet(x_point, y_point, cur_yaw, map_.x, map_.y);
    next_s_vals_.push_back(cur_sd[0]);
    next_d_vals_.push_back(d);
    if (DEBUG) {
      std::cout << "Added new trajectory step " << (t+1) << " with x: " << x_point << " y: " << y_point <<
          " keep_speed: 1" << std::endl;
    }
  }
}


/**
 * Method to apply the pre-computed trajectory or generate new one if one is not computed
 * @param p The latest position data planned
 * @param forward_traj The pre-computed forward trajectory
 * @param cur_traj_idx The current position on the forward_traj trah=jectory
 */
void FSM::applyTrajectory(PositionData & p, vector<TrajNode> & forward_traj, size_t & cur_traj_idx)
{
  updateLane(p.d);
  createXYFromSpline(p.d, forward_traj, cur_traj_idx);
}

/**
 * Method to clear processed data from the FSM
 */
void FSM::clearProcessed()
{
  if (first_time_) return; // No processed to clear during the first time we run
  int processed = TRAJ_STEPS - cxt_->path_x.size();
  assert(processed > 0 and processed <= TRAJ_STEPS);

  if (DEBUG) {
    for (int i = 0; i < processed; ++i) {
      std::cout << (i+1) << " Processed x: " << next_x_vals_[i] << " y: " << next_y_vals_[i] << std::endl;
    }
  }
  next_x_vals_.erase(next_x_vals_.begin(), next_x_vals_.begin() + processed);
  next_y_vals_.erase(next_y_vals_.begin(), next_y_vals_.begin() + processed);
  next_s_vals_.erase(next_s_vals_.begin(), next_s_vals_.begin() + processed);
  next_d_vals_.erase(next_d_vals_.begin(), next_d_vals_.begin() + processed);
}

/**
 * Method to get the last data that was processed by the FSM (x,y,s,d)
 * @return The last processed PositionData
 */
PositionData FSM::getLastProcessed() const
{
  PositionData p;
  if (first_time_) {
    p.x = cxt_->x; p.y = cxt_->y; p.s = cxt_->s;
    p.d = STARTING_LANE;
  } else {
    int processed = TRAJ_STEPS - cxt_->path_x.size();
    assert(processed > 0 and processed <= TRAJ_STEPS);

    p.x = next_x_vals_[processed -1];
    p.y = next_y_vals_[processed -1];
    p.s = next_s_vals_[processed -1];
    p.d = next_d_vals_[processed -1];
  }
  return p;
}

/**
 * Method to get the most recent data that on the processing queue of the FSM (x,y,s,d)
 * @return The most recently planned position data
 */
PositionData FSM::getMostRecentNotProcessed() const
{
  PositionData p;
  if (first_time_) {
    if (DEBUG) std::cout << "getting position data for first time\n";
    p.x = cxt_->x; p.y = cxt_->y; p.s = cxt_->s;
    p.d = STARTING_LANE;
  } else {
    p.x = next_x_vals_.back();
    p.y = next_y_vals_.back();
    p.s = next_s_vals_.back();
    p.d = next_d_vals_.back();
  }
  return p;
}

/**
 * Method to handle first time use. Will just accelerate to MAX_SPEED * MAX_SPEED_MAR
 */
void FSM::handleFirstTime()
{
  assert(state_ == nullptr);
  state_ = new MatchSpeed(*this, MAX_SPEED * MAX_SPEED_MAR, nullptr);
  cur_policy_ = MATCH_SPEED;
  first_time_ = false;
}

/**
 * Main update function for the FSM based on new State
 * @param cxt Context containing the latest data for the state of the car
 */
void FSM::update(Context & cxt) 
{
  cxt_ = &cxt;
  if (first_time_) {
    if (DEBUG) std::cout << "Running for the first time\n";
    handleFirstTime();
    return;
  }

  int sz = cxt.path_x.size();
  // How many points left in the previous path?
  if (sz == TRAJ_STEPS) {
    if (DEBUG) std::cout << "Returning same lists since nothing was processed\n";
    return; // nothing processed, the lists remain the same
  }

  PositionData p = getLastProcessed();

  if (DEBUG) {
    static int stepCnt = 0;
    std::cout << "================================================================================\n";
    std::cout << "[Step " << stepCnt << "] Current given state  (s,d,x,y,yaw,v): (" <<
        cxt_->s << ", " << cxt_->d << ", " << cxt_->x << ", " << cxt_->y << ", " << deg2rad(cxt_->yaw) <<
        ", " << cxt_->speed * MPH2MPS << ")\n";
    std::cout << "[Step " << stepCnt << "] Current internal state  (s,d,x,y,yaw): (" <<
        p.s << ", " << p.d << ", " << p.x << ", " << p.y << ", " << deg2rad(cxt_->yaw) << ")\n";
    std::cout << "State target speed: : " << state_->getTargetSpeed() << std::endl;
    stepCnt += (TRAJ_STEPS - sz);
    std::cout << "FSM update. Server processed " << (TRAJ_STEPS - sz) << " points" << std::endl;
  }

  Policy policy = state_->nextPolicy();

  // Are we already in the process of a maneuver? Then go straight to update
  if (policy == KEEP_LANE) {
    // Check sensor fusion to see if we need to change lane
    SensorData sd = checkLanes(p);
    // Check if our policy is changed and we need to change lane
    if (sd.policy == CHANGE_LANE) {
      delete state_;
      state_ = new ChangeLane(*this, sd.v, &sd);
      return;
    // Otherwise, check if our speed target has changed
    } else if (state_->getTargetSpeed() - sd.v > 0.1 or state_->getTargetSpeed() - sd.v < -0.1) {
      delete state_;
      state_ = new MatchSpeed(*this, sd.v, &sd);
      return;
    // If we just now changed to KEEP_LANE, change the state
    } else if (policy == KEEP_LANE and cur_policy_ != policy) {
      delete state_;
      state_ = new KeepLane(*this, sd.v);
    }
  }

  cur_policy_ = policy;

  // Update current state
  state_->update();
  return;
}
