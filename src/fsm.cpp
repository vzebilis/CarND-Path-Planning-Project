#include <iostream>
#include <cassert>
#include "fsm.h"
#include "helpers.h"

using nlohmann::json;
using std::string;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;


static double getBestStep(vector<double> & coeffs, double T) 
{
  double T2 = T * T;
  double T3 = T2 * T;
  double T4 = T3 * T;
  double T5 = T4 * T;
  return coeffs[0] + coeffs[1]*T + coeffs[2]*T2 + coeffs[3]*T3 + coeffs[4]*T4 + coeffs[5]*T5;
}

// Function to compute the move on s with lowest cost in jerk
// Only using the s, speed, and acc of given states, over the given time range t
static vector<double> computePolyCoefficients(TrajData & traj_data)
{
  double es   = traj_data.final_s; // End s, speed and acc
  double esd  = traj_data.final_speed;
  double esdd = traj_data.final_acc;
  double ss   = traj_data.init_s; // Start s, speed and acc
  double ssd  = traj_data.init_speed;
  double ssdd = traj_data.init_acc;
  double t    = traj_data.time;
  double t2   = t * t;
  double t3   = t2 * t;
  double t4   = t3 * t;
  double t5   = t4 * t;
  MatrixXd TMat(3, 3);
  TMat <<    t3,    t4,    t5,
           3*t2,  4*t3,  5*t4,
            6*t, 12*t2, 20*t3;
  MatrixXd TMat_i = TMat.inverse();

  VectorXd SVec(3);
  SVec << es - (ss + ssd * t + ssdd * t2 * 0.5),
          esd - (ssd + ssdd * t),
          esdd - ssdd;

  VectorXd a345 = TMat_i * SVec;  // result is a vector of 3 elements
  
  return {ss, ssd, ssdd * 0.5, a345[0], a345[1], a345[2]};
}


// Transform from Frenet s coordinate to Cartesian x,y using Splines.
// d_offset is the lane offset (2: left lane, 6: middle lane, 10: right lane)
vector<double> FSM::getXYfromSpline(double s, double d_offset)
{
  int prev_wp = -1;

  while (s > map_.s[prev_wp+1] && (prev_wp < (int)(map_.s.size()-1))) {
    ++prev_wp;
  }

  int wp2 = (prev_wp+1)%map_.x.size();

  double prev_x = spd_.x(map_.s[prev_wp]);
  double prev_y = spd_.y(map_.s[prev_wp]);
  double x2 = spd_.x(map_.s[wp2]);
  double y2 = spd_.y(map_.s[wp2]);
  double heading = atan2(y2 - prev_y, x2 - prev_x);

  double seg_x  = spd_.x(s);
  double seg_y  = spd_.y(s);
  //double seg_dx = spd_.sdx(s);
  //double seg_dy = spd_.sdy(s);
  //double d = sqrt(pow(seg_dx, 2) + pow(seg_dy, 2)) * d_offset; // Chosen lane
  double d = d_offset; // Chosen lane

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};
}

FSM::FSM(const Map & map) : 
  state_(nullptr), first_time_(true), cxt_(nullptr), prev_pol_(KEEP_LANE), map_(map),
  prev_trg_speed_(0), delta_speed_to_zero_acc_(0)
{
  static_assert(STARTING_LANE == 2 or STARTING_LANE == 6 or STARTING_LANE == 10);
  // Set current lane based on STARTING_LANE constexpr
  if (STARTING_LANE == 2) lane_ = LEFT;
  else if (STARTING_LANE == 6) lane_ = CENTER;
  else lane_ = RIGHT;

  // Generate the full track splines
  generateFullSplines();

  // Compute some const information
  int shed_acc_steps = ceil((MAX_ACC / MAX_JERK) / TIME_RES);
  double a_shed = 0;
  for (int i = 0; i < shed_acc_steps; ++i) {
    delta_speed_to_zero_acc_ += a_shed * TIME_RES;
    a_shed += MAX_JERK * TIME_RES;
  }
}

// Generate the full splines based on map data
// Will be used to compute trajectories
void FSM::generateFullSplines()
{
  spd_.x.set_points(map_.s, map_.x);
  spd_.y.set_points(map_.s, map_.y);
  spd_.dx.set_points(map_.s, map_.dx);
  spd_.dy.set_points(map_.s, map_.dy);
}

// Method to compute and return the minimum amount of road length and time
// it would take to acceleratate to target speed with 0 final acceleration
TrajData FSM::computeMatchTargetSpeed(double init_s, double init_d, double init_speed, double init_acc, 
                                      double final_d, double trg_speed) 
{
  assert(init_s     >= 0);
  assert(init_speed >= 0);
  assert(trg_speed  <= MAX_SPEED);

  if (DEBUG) std::cout << "Try to match target speed: " << trg_speed << std::endl;

  TrajData ret; 
  ret.final_speed         = init_speed; // default return value
  ret.init_d              = init_d;     
  ret.final_d             = final_d;     

  double delta_v = trg_speed - init_speed;
  // If we are already at the target speed, with zero acceleration, then nothing to do
  if (delta_v > 0 and delta_v < 0.1 and isZero(init_acc)) return ret;


  // How many timesteps do we need at least to shed our initial acceleration?
  int shed_acc_steps      = ceil((fabs(init_acc) / MAX_JERK) / TIME_RES);
  double delta_speed_shed = 0;
  double a_shed           = 0;
  double a_jerk           = (init_acc > 0)? -MAX_JERK: MAX_JERK;
  // How much speed are we going to add while shedding initial acceleration?
  for (int i = 0; i < shed_acc_steps; ++i) {
    a_shed               += a_jerk * TIME_RES;
    delta_speed_shed     += a_shed * TIME_RES;
  }
  const bool accel = delta_v + delta_speed_shed >= 0;

  // The flip point: speed at which to switch from accelerating to decelerating or vice versa
  double flip_point = init_speed + (delta_v + delta_speed_shed) / 2;
  double T = 0;
  double S = init_s;
  double cur_speed = init_speed;
  double cur_acc = init_acc;
  while (((accel)? trg_speed > cur_speed: cur_speed > trg_speed)) {
    if (DEBUG) {
      std::cout << "Cur speed: " << cur_speed << " Cur acc: " << cur_acc << " Cur S: " << S <<
                   " Speed flip point: " << flip_point << std::endl;
    }
    T           += TIME_RES;
    S           += cur_speed * TIME_RES;
    cur_speed   += cur_acc * TIME_RES;
    bool flipped = (accel)? cur_speed > flip_point: cur_speed < flip_point;
    double jerk  = (accel)? MAX_JERK: -MAX_JERK;
    jerk         = (flipped)? -jerk: jerk;
    cur_acc     += jerk * TIME_RES;
    if (fabs(cur_acc) > MAX_ACC) {
      cur_acc = (cur_acc > 0)? MAX_ACC: -MAX_ACC;
      // Update flip point as well
      flip_point = trg_speed + ((cur_acc > 0)? -1: 1) * delta_speed_to_zero_acc_;
    }
    if (flipped and isZero(cur_acc)) break;
  }

  ret.init_s      = init_s;
  ret.init_speed  = init_speed;
  ret.init_acc    = init_acc;
  ret.final_s     = S;
  ret.final_speed = trg_speed;
  ret.final_acc   = 0;
  ret.time        = T * 1.1; // Add 10% margin

  return ret;
}

SplineData * FSM::getShortSplines(TrajData & td)
{
  SplineData * spd =  new SplineData;
  int msz = map_.s.size();
  assert(td.time > 0);
  assert(map_.s.size() > 0);

  // Get next WP from final position
  int next_wp = msz;
  while (td.final_s < map_.s[next_wp-1] and (next_wp > 0)) --next_wp;
  int prev_wp = -1;
  while (td.init_s > map_.s[prev_wp+1] and (prev_wp < msz-1)) ++prev_wp;
  // Sanity checks
  assert(prev_wp > 0);
  assert(next_wp < msz);

  auto next_s     = map_.s[next_wp];
  double next_t   = next_s - fmod(td.final_s, MAX_S);
  if (next_t < 0) next_t += MAX_S;
  next_t /= td.final_speed? td.final_speed: 0.01;
  next_t += td.time;

//  auto nnext_wp = next_wp +1;
//  nnext_wp = nnext_wp % msz;
//  auto nnext_s    = map_.s[nnext_wp];
//  double nnext_t  = nnext_s - next_s;
//  if (nnext_t < 0) nnext_t += MAX_S;
//  nnext_t /= td.final_speed? td.final_speed: 0.01;
//  nnext_t += next_t;
  // Use beginning of current trajectory data as previous waypoint
//  auto prev_s     = next_s_vals_.front();
//  double prev_t   = fmod(td.init_s, MAX_S) - fmod(prev_s, MAX_S);
//  if (prev_t < 0) prev_t += MAX_S;
//  prev_t /= td.init_speed? td.init_speed: 0.01;

  auto prev_s     = map_.s[prev_wp];
  double prev_t   = fmod(td.init_s, MAX_S) - prev_s;
  if (prev_t < 0) prev_t += MAX_S;
  prev_t /= td.init_speed? td.init_speed: 0.01;

//  auto pprev_wp = prev_wp -1;
//  pprev_wp = pprev_wp % msz;
//  auto pprev_s     = map_.s[pprev_wp];
//  double pprev_t   = prev_s - pprev_s;
//  if (pprev_t < 0) pprev_t += MAX_S;
//  pprev_t /= td.init_speed? td.init_speed: 0.01;
//  pprev_t += prev_t;

  if (DEBUG) std::cout << "GetShortSplines. time: " << td.time << " init_s: " << td.init_s << " final_s: " <<
      td.final_s << " next_s: " << next_s << " next_t: " << next_t << " prev_s: " << prev_s << " prev_t: " <<
      -prev_t << std::endl;

  auto prev_xy  = getXY(prev_s, td.init_d, map_.s, map_.x, map_.y);
  auto init_xy  = getXY(td.init_s, td.init_d, map_.s, map_.x, map_.y);
  auto final_xy = getXY(td.final_s, td.final_d, map_.s, map_.x, map_.y);
  auto next_xy  = getXY(next_s, td.final_d, map_.s, map_.x, map_.y);

  vector<double> times = {-prev_t, 0, td.time, next_t};
  spd->x.set_points(times, {prev_xy[0], init_xy[0], final_xy[0], next_xy[0]});
  spd->y.set_points(times, {prev_xy[1], init_xy[1], final_xy[1], next_xy[1]});
  spd->s.set_points(times, {prev_s, td.init_s, td.final_s, next_s});
  spd->d.set_points(times, {td.init_d, td.init_d, td.final_d, td.final_d});

  if (false and DEBUG) {
    std::cout << "GetShortSplines. Going through all the steps in the spline s:" <<std::endl;
    double time = TIME_RES;
    while (time <= td.time) {
      std::cout << "Time: " << time << " S: " << spd->s(time) << std::endl;
      time += TIME_RES;
    }
  }
  return spd;
}

//
//
//
void ChangeLane::computeTrajectory(TrajData & td)
{
  forward_traj_.clear();
  for (size_t i = 0; i < fsm_.next_x_vals_.size(); ++i) forward_traj_.emplace_back();
  cur_traj_idx_  = 0;

  auto spd = fsm_.getShortSplines(td);

  // Generate trajectory until we reach the required time
  double time    = TIME_RES;
  while (time <= td.time) {
    forward_traj_.emplace_back();
    auto & ft    = forward_traj_.back();
    ft.x         = spd->x(time);
    ft.y         = spd->y(time);
    ft.s         = spd->s(time);
    ft.d         = spd->d(time);
    time        += TIME_RES;
  }
  delete spd;
  if (DEBUG) {
    std::cout << "Computed new trajectory with " << forward_traj_.size() << " steps to target speed " <<
            trg_speed_ << " in time " << td.time << " :\n";
    std::cout << " init_s: " << td.init_s << " init_d: " << td.init_d << " init_speed: " << td.init_speed <<
            " init_acc: " << td.init_acc << " final_s: " << td.final_s << " final_d: " << td.final_d <<
            " final_speed: " << td.final_speed << " final_acc: " << td.final_acc << std::endl;
  }
}

//
// Generic state method to compute the necessary trajectory based on
// the supplied TrajData info
//
void State::computeTrajectory(TrajData & td)
{
  assert(td.time);
  size_t sz = cur_traj_idx_ + fsm_.cxt_->path_x.size();
  if (sz > forward_traj_.size()) forward_traj_.clear();
  else forward_traj_.erase(forward_traj_.begin() + sz, forward_traj_.end());
  //forward_traj_.clear();
  //cur_traj_idx_  = 0;
  auto coeffs    = computePolyCoefficients(td);
  // Generate trajectory until we reach the required time
  double time    = TIME_RES;
  // Calculate step for D if there is a significant difference
  double diff_d  = td.final_d - td.init_d;
  double step_d  = (fabs(diff_d) < 3)? 0: diff_d / (td.time / TIME_RES);
  if (DEBUG) std::cout << "computeTrajectory. diff_d: " << diff_d << " step_d: " << step_d << " time: " <<  td.time << std::endl;
  double new_d   = td.init_d + step_d;
  while (time <= td.time) {
    forward_traj_.emplace_back();
    auto & ft    = forward_traj_.back();
    double new_s = getBestStep(coeffs, time);
    new_s        = fmod(new_s, MAX_S);
    auto vecXY   = fsm_.getXYfromSpline(new_s, new_d);
    ft.x         = vecXY[0];
    ft.y         = vecXY[1];
    ft.s         = new_s;
    ft.d         = new_d;
    time        += TIME_RES;
    new_d       += step_d;
  }

  if (DEBUG) {
    std::cout << "State. Computed new trajectory with " << forward_traj_.size() << " steps to target speed " <<
                 trg_speed_ << " in time " << td.time << " :\n";
    std::cout << " init_s: " << td.init_s << " init_d: " << td.init_d << " init_speed: " << td.init_speed <<
                 " init_acc: " << td.init_acc << " final_s: " << td.final_s << " final_d: " << td.final_d <<
                 " final_speed: " << td.final_speed << " final_acc: " << td.final_acc << std::endl;
  }
}


//
// Function to return the closest car ahead or behind, in the lane defined by d_offset
//
static SensorData getSensorFusion(const PositionData & p, Context & cxt, bool ahead, double d_offset)
{
  static constexpr double HALF_LANE_WIDTH = LANE_WIDTH / 2;

  // If we are looking at the side lanes, look a bit further ahead (not behind!)
  double factor = 1.0;
  if (d_offset > HALF_LANE_WIDTH or d_offset < -HALF_LANE_WIDTH) factor = 1.5;

  // Get data for car closest to us in s, going forward
  int min_idx  = -1;
  double min_diff_s = MAX_S;
  for (int i = 0; i < cxt.sensor_fusion.size(); ++i) {
    // Only interested about cars in the same lane for the given offset
    double cur_d = cxt.sensor_fusion[i][6]; // get the d value
    if (cur_d < p.d + d_offset - HALF_LANE_WIDTH or
        cur_d > p.d + d_offset + HALF_LANE_WIDTH) continue;

    double st_s   = p.s;
    double cur_s  = cxt.sensor_fusion[i][5]; // get the s value

    // Handle wraparound
    if (ahead) {
      if (cur_s < st_s) cur_s += MAX_S;
    } else {
      if (cur_s > st_s) st_s  += MAX_S; 
    }
    
    // Only interested about cars in a reasonable distance
    double diff_s = fabs(cur_s - st_s);
    double margin = SENSOR_MAX_DIST;
    if (cur_s - st_s > 0) margin *= factor;
    if (diff_s > margin) continue;

    if (diff_s < min_diff_s) {
      min_diff_s  = diff_s;
      min_idx     = i;
    }
  }
  SensorData sd; 
  sd.diff_s = min_diff_s;
  // Are there any cars detectable in front of us?
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

static std::pair<double, double> getTimeToIntercept(double dist, double speed, double trg_speed) 
{
  assert(speed > trg_speed);
  double time = 0;
  double trav = 0;
  while (dist > 0) {
    trav += speed * TIME_RES;
    dist -= (speed - trg_speed) * TIME_RES;
    time += TIME_RES;
  }
  return {time, trav};
}

// Function that decides whether we can change lane, based on cars ahead and behind
bool FSM::canChangeLane(const PositionData & p, SensorData & ahead, SensorData & behind)
{
  double time_steps_before_change = cxt_->path_x.size() * TIME_RES;
  double time_to_change = LANE_CHANGE_FCT / p.v; // time to change lane with current speed
  // First check the car behind
  {
    if (behind.diff_s < LANE_CHANGE_MAR) return false;
    double dv = behind.v - p.v;
    double future_diff_s = behind.diff_s - LANE_CHANGE_MAR - time_steps_before_change * dv;
    if (future_diff_s < dv * time_to_change) return false;
  }
  // Now check the car ahead
  {
    if (ahead.diff_s < LANE_CHANGE_MAR) return false;
    double dv = p.v - ahead.v;
    double future_diff_s = ahead.diff_s - LANE_CHANGE_MAR - time_steps_before_change * dv;
    if (future_diff_s < dv * time_to_change) return false;
  }
  return true;
}
//
// Function to check which of the available lanes is best to use
//
SensorData FSM::checkLanes(const PositionData & p)
{
  if (DEBUG) std::cout << "checkLanes: current lane: " << lane_ << std::endl;
  // Check same lane
  SensorData best = getSensorFusion(p, *cxt_, true, 0); best.d = p.d; best.pol = KEEP_LANE;

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
    auto right_ahead = getSensorFusion(p, *cxt_, true,  4);
    auto right_backw = getSensorFusion(p, *cxt_, false, 4);
    // Check if there is enough space
    if (canChangeLane(p, right_ahead, right_backw)) {
      // If there are no cars present ahead right or if their speed is higher than the best, chose right
      if (right_ahead.s == -1 or best.v < 0.95 * right_ahead.v) {
        best = right_ahead;
        best.pol = CHANGE_LANE;
        best.d = (lane_ == LEFT)? 6: 10;
        if (best.s == -1)  {
          if (DEBUG) std::cout << "Moving to right lane, since there is no car there\n";
          best.v = prev_trg_speed_; //MAX_SPEED * MAX_SPEED_MAR;
          return best;
        }
        if (DEBUG) std::cout << "Car in our right lane: " << best.id << " s: " << best.s << " d: " << best.d << 
                     " v: " << best.v << std::endl;
      }
    }
  }

  // Check lane to our left
  if (lane_ == RIGHT or lane_ == CENTER) {
    auto left_ahead = getSensorFusion(p, *cxt_, true,  -4);
    auto left_backw = getSensorFusion(p, *cxt_, false, -4);
    // Check if there is enough space
    if (canChangeLane(p, left_ahead, left_backw)) {
      // If there are no cars present ahead left or if their speed is higher than the best, chose left
      if (left_ahead.s == -1 or best.v < 0.95 * left_ahead.v) {
        best = left_ahead;
        best.pol = CHANGE_LANE;
        best.d = (lane_ == RIGHT)? 6: 2;
        if (best.s == -1)  {
          if (DEBUG) std::cout << "Moving to left lane, since there is no car there\n";
          best.v = prev_trg_speed_; //MAX_SPEED * MAX_SPEED_MAR;
          return best;
        }
        if (DEBUG) std::cout << "Car in our left lane: " << best.id << " s: " << best.s << " d: " << best.d << 
                     " v: " << best.v << std::endl;
      }
    }
  }

  if (DEBUG) std::cout << "Applying policy: " << best.pol << std::endl;
  return best;
}

//
// Method to update the internal lane state
//
void FSM::updateLane(double d)
{
  assert(d >= 0 and d <= 12);
  if (d >= 0 and d < 4) lane_ = LEFT;
  else if (d >= 4 and d < 8) lane_ = CENTER;
  else if (d >= 8 and d < 12) lane_ = RIGHT;
}

//
// Change Lane STATE
//
ChangeLane::ChangeLane(FSM & fsm, double trg_speed, SensorData * sd) : State(fsm, trg_speed)
{
  if (DEBUG) std::cout << "ChangeLane: creating new State\n";
  //p_ = fsm_.getLastProcessed();
  //fsm_.clearAll();
  p_ = fsm_.getMostRecentNotProcessed();
  fsm_.clearProcessed();
  TrajData td = computeTargetPos(p_, sd);
  fsm_.prev_trg_speed_ = trg_speed_;
  computeTrajectory(td); // Base class method
  // Apply the new trajectory by creating the new next_x_vals and next_y_vals
  fsm_.applyTrajectory(p_, forward_traj_, cur_traj_idx_);
}

//
// Compute target position, for a lane change
//
TrajData ChangeLane::computeTargetPos(PositionData & p, SensorData * sd)
{
  TrajData td;
  td.init_s      = p.s;
  td.init_d      = p.d;
  td.init_speed  = p.v;
  td.init_acc    = p.a;

  // If changing lanes complete lane change fast
  td.time        = LANE_CHANGE_FCT / p.v; // time to switch lanes
  td.final_s     = p.s + p.v * td.time * 1.1; // distance travelled
  td.final_d     = sd->d;
  td.final_speed = p.v;
  td.final_acc   = 0;

  td.mid_s       = (td.init_s + td.final_s) / 2;
  td.mid_d       = (td.init_d + td.final_d) / 2;

  if (DEBUG) std::cout << "Changing lanes from d: " << td.init_d << " to d: " << td.final_d << " in time: " << td.time <<
          " and s: " << (p.v * td.time) << std::endl;
  return td;
}

//
// Update method to be called for every update of the FSM
//
void ChangeLane::update()
{
  p_ = fsm_.getMostRecentNotProcessed();
  // Clear old processed data on the FSM
  fsm_.clearProcessed();
  fsm_.incrementByProcessed(cur_traj_idx_);
  // Apply the new trajectory by creating the new next_x_vals and next_y_vals
  fsm_.applyTrajectory(p_, forward_traj_, cur_traj_idx_);
}

//
// Method determining next policy
//
Policy ChangeLane::nextPolicy()
{
  if (cur_traj_idx_ >= forward_traj_.size()) {
    if (DEBUG) std::cout << "ChangeLane: next policy is KEEP_LANE\n";
    return KEEP_LANE;
  }
  if (DEBUG) std::cout << "ChangeLane: next policy remains CHANGE_LANE\n";
  return CHANGE_LANE;
}


//
// Keep Lane (and Speed) update
//
void KeepLane::update()
{
  p_ = fsm_.getMostRecentNotProcessed();
  // Clear old processed data on the FSM
  fsm_.clearProcessed();
  fsm_.applyTrajectory(p_, forward_traj_, cur_traj_idx_);
}

//
// Match Speed STATE
//
MatchSpeed::MatchSpeed(FSM & fsm, double trg_speed, SensorData * sd) : State(fsm, trg_speed)
{
  if (DEBUG) std::cout << "MatchSpeed: creating new State\n";
  p_ = fsm_.getLastProcessed();
  //p_ = fsm_.getMostRecentNotProcessed();
  fsm_.clearAll();
  TrajData td = computeTargetPos(p_, sd);
  fsm_.prev_trg_speed_ = trg_speed_;
  computeTrajectory(td); // Base class method
  // Apply the new trajectory by creating the new next_x_vals and next_y_vals
  fsm_.applyTrajectory(p_, forward_traj_, cur_traj_idx_);
}

//
// Update method to be called for every update of the FSM
//
void MatchSpeed::update()
{
  p_ = fsm_.getMostRecentNotProcessed();
  // Clear old processed data on the FSM
  fsm_.clearProcessed();
  fsm_.incrementByProcessed(cur_traj_idx_);
  // Apply the new trajectory by creating the new next_x_vals and next_y_vals
  fsm_.applyTrajectory(p_, forward_traj_, cur_traj_idx_);
}


//
// Compute target position, for matching observed speed
//
TrajData MatchSpeed::computeTargetPos(PositionData & p, SensorData * sd)
{
  // If no car spotted, go to max speed
  if (!sd or sd->s == -1) {
    if (DEBUG) std::cout << "Setting target speed to MAX\n";
    return fsm_.computeMatchTargetSpeed(p.s, p.d, p.v, p.a, p.d, MAX_SPEED * MAX_SPEED_MAR);
  }

  // Compute target speed
  double trg_speed = sd->v;
  trg_speed        = fmin(trg_speed, MAX_SPEED * MAX_SPEED_MAR);

  TrajData td      = fsm_.computeMatchTargetSpeed(p.s, p.d, p.v, p.a, sd->d, trg_speed);
  // Accelerate to match the preceding car's speed
  // Don't waste time
  if (trg_speed >= p.v) return td;

  // Otherwise, decelerate more smoothly if we have enough space/time to do it
  // In order to so so, we calculate the future position of the observed car
  // in the time it will take us to match the target speed.
  // This give us an extra margin to decelerate more smoothly
  double car_s = sd->s - FRONT_CAR_MAR;
  if (car_s < p.s) car_s += MAX_S;
  double car_future_s = car_s + td.time * trg_speed;

  // If there is no time to stop before a collision, violate constraints to do so properly
  if (td.final_s > car_future_s) {
    td.final_s = car_future_s;
  // Otherwise just decelerate as smoothly as possible
  } else {
    td.time += (car_future_s - td.final_s) / p.v;
    td.final_s = car_future_s;
  }

  // TODO: we should use the middle time step value of s
  td.mid_s = (td.init_s + td.final_s) / 2;
  td.mid_d = p.d;
  //auto time_trav = getTimeToIntercept(car_s - init_s, init_speed, trg_speed);
  //if (DEBUG) std::cout << "Time to intercept car: " << time_trav.first << " and s: " << time_trav.second << std::endl;
  //td.final_s     = init_s + time_trav.second;
  //td.final_speed = trg_speed;
  //td.time        = time_trav.first;
  return td;
}

//
// Method determining next policy
//
Policy MatchSpeed::nextPolicy()
{
  if (cur_traj_idx_ >= forward_traj_.size()) {
    if (DEBUG) std::cout << "MatchSpeed: next policy is KEEP_LANE\n";
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
void FSM::incrementByProcessed(size_t & val) const
{
  if (first_time_) return;
  size_t processed = TRAJ_STEPS - cxt_->path_x.size();
  assert(processed > 0 and processed <= TRAJ_STEPS);
  val += processed;
}

static void printPositionData(std::string prefix, PositionData & p)
{
  std::cout << prefix.c_str() << " s: " << p.s << " d: " << p.d << " v: " << p.v << " a: " << p.a << std::endl;
}

static TrajData initTrajData(size_t sz, double init_s, double init_speed, double init_d)
{
  TrajData td;
  td.time = (TRAJ_STEPS - sz) * TIME_RES;
  td.init_s = td.final_s = init_s;
  td.final_s += td.init_speed * td.time;
  td.init_speed = td.final_speed = init_speed;
  td.init_acc = td.final_acc = 0;
  td.init_d = td.final_d = td.mid_d = init_d;
  td.mid_s = (td.init_s + td.final_s) / 2;
  return td;
}

//
// Method to apply the pre-computed trajectory
//
void FSM::applyTrajectory(PositionData & p, vector<TrajNode> & forward_traj, size_t cur_traj_idx)
{
  double prev_ns   = p.s;
  double prev_d    = p.d;
  double cur_speed = p.v;
  double cur_acc   = p.a;

  double prev_speed = cur_speed;

  updateLane(prev_d);

  SplineData * spd = nullptr;
  bool keep_speed = (cur_traj_idx >= forward_traj.size());

  size_t i = next_x_vals_.size();

//  if (keep_speed) {
//    TrajData td = initTrajData(i, prev_ns, prev_speed, prev_d);
//    spd = getShortSplines(td);
//  }

  if (DEBUG) printPositionData("ApplyTrajectory", p);

  double time = 0;
  // Build the rest of the steps
  for (; i < TRAJ_STEPS; ++i) {
    // If we completed our computed trajectory, just change policy to keep current speed
    if (!keep_speed and cur_traj_idx + i >= forward_traj.size()) {
      //TrajData td = initTrajData(i, prev_ns, prev_speed, prev_d);
      //spd = getShortSplines(td);
      keep_speed = true;
      time = 0;
      if (DEBUG) std::cout << "Finished generated trajectory. Now keeping speed (TrajSize: " << forward_traj.size() << ")\n";
    }
    time += TIME_RES;
    double nx, ny, ns, nd;
    if (keep_speed) {
//      nx = spd->x(time);
//      ny = spd->y(time);
//      ns = spd->s(time);
//      nd = prev_d;
      ns = prev_ns + prev_speed * TIME_RES;  // TODO: perhaps use prev_speed?
      ns = fmod(ns, MAX_S);
      nd = prev_d;
      auto vecXY = getXYfromSpline(ns, nd);
      nx = vecXY[0];
      ny = vecXY[1];
    } else {
      nx = forward_traj[cur_traj_idx + i].x;
      ny = forward_traj[cur_traj_idx + i].y;
      ns = forward_traj[cur_traj_idx + i].s;
      nd = forward_traj[cur_traj_idx + i].d;
    }
    if (DEBUG) {
      std::cout << "Added new trajectory step " << (i+1) << " with s: " << ns << " d: " << nd << " speed: " << cur_speed <<
          " acc: " << cur_acc << " keep_speed: " << keep_speed << std::endl;
    }
    next_x_vals_.push_back(nx);
    next_y_vals_.push_back(ny);
    next_s_vals_.push_back(ns);
    next_d_vals_.push_back(nd);
    // Keep track of speed and acceleration
    if (ns < prev_ns) ns += MAX_S; // Wraparound case
    cur_speed  = (ns - prev_ns) / TIME_RES;
    cur_speed  = fmin(cur_speed, MAX_SPEED * MAX_SPEED_MAR);
    cur_acc    = (cur_speed - prev_speed) / TIME_RES;
    cur_acc    = fmin(cur_acc, MAX_ACC);
    prev_ns    = fmod(ns, MAX_S);
    prev_d     = nd;
    prev_speed = cur_speed;
    prev_acc_.push_back(cur_acc);
    prev_speed_.push_back(cur_speed);
  }
  if (spd) delete spd;
}

//
// Method to clear processed data from the FSM
//
void FSM::clearProcessed()
{
  if (first_time_) return; // No processed to clear during the first time we run
  int processed = TRAJ_STEPS - cxt_->path_x.size();
  assert(processed > 0 and processed <= TRAJ_STEPS);

  next_x_vals_.erase(next_x_vals_.begin(), next_x_vals_.begin() + processed);
  next_y_vals_.erase(next_y_vals_.begin(), next_y_vals_.begin() + processed);
  next_s_vals_.erase(next_s_vals_.begin(), next_s_vals_.begin() + processed);
  next_d_vals_.erase(next_d_vals_.begin(), next_d_vals_.begin() + processed);
  prev_acc_.erase(prev_acc_.begin(), prev_acc_.begin() + processed);
  prev_speed_.erase(prev_speed_.begin(), prev_speed_.begin() + processed);
}

//
// Method to get the last data that was processed by the FSM (x,y,s,d,v,a)
//
PositionData FSM::getLastProcessed() const
{
  PositionData p;
  if (first_time_) {
    p.x = cxt_->x; p.y = cxt_->y; p.s = cxt_->s;
    p.d = STARTING_LANE; p.v = 0; p.a = 0;
  } else {
    int processed = TRAJ_STEPS - cxt_->path_x.size();
    assert(processed > 0 and processed <= TRAJ_STEPS);

    p.x = next_x_vals_[processed -1];
    p.y = next_y_vals_[processed -1];
    p.s = next_s_vals_[processed -1];
    p.d = next_d_vals_[processed -1];
    p.v = prev_speed_[processed -1];
    p.a = prev_acc_[processed -1];
  }
  return p;
}

//
// Method to get the most recent data that on the processing queue of the FSM (x,y,s,d,v,a)
//
PositionData FSM::getMostRecentNotProcessed() const
{
  PositionData p;
  if (first_time_) {
    if (DEBUG) std::cout << "getting position data for first time\n";
    p.x = cxt_->x; p.y = cxt_->y; p.s = cxt_->s;
    p.d = STARTING_LANE; p.v = 0; p.a = 0;
  } else {
    p.x = next_x_vals_.back();
    p.y = next_y_vals_.back();
    p.s = next_s_vals_.back();
    p.d = next_d_vals_.back();
    p.v = prev_speed_.back();
    p.a = prev_acc_.back();
  }
  return p;
}

//
// Method to handle first time use. Will just accelerate to MAX_SPEED
//
void FSM::handleFirstTime()
{
  assert(state_ == nullptr);
  state_ = new MatchSpeed(*this, MAX_SPEED * MAX_SPEED_MAR, nullptr);
  first_time_ = false;
}

//
// Main update function for the FSM based on new State
//
void FSM::update(Context & cxt) 
{
  cxt_ = &cxt;
  if (first_time_) {
    if (DEBUG) std::cout << "Running for the first time\n";
    handleFirstTime();
    prev_pol_ = MATCH_SPEED;
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
    std::cout << "[Step " << stepCnt << "] Current internal state  (s,d,x,y,v,a): (" <<
        p.s << ", " << p.d << ", " << p.x << ", " << p.y << ", " << p.v << ", " << p.a << ")\n";
    std::cout << "State target speed: : " << state_->getTargetSpeed() << std::endl;
    stepCnt += (TRAJ_STEPS - sz);
    std::cout << "FSM update. Server processed " << (TRAJ_STEPS - sz) << " points" << std::endl;
  }

  Policy pol = state_->nextPolicy();

  // Are we already in the process of a maneuver? Then go straight to update
  if (pol == KEEP_LANE) {
    // Check sensor fusion to see if we need to change lane
    SensorData sd = checkLanes(p);
    // Check if our policy is changed and we need to change lane
    if (sd.pol == CHANGE_LANE) {
      pol = CHANGE_LANE;
      delete state_;
      state_ = new ChangeLane(*this, sd.v, &sd);
      prev_pol_ = pol;
      return;
    // Otherwise, check if our speed target has changed
    } else if (state_->getTargetSpeed() - sd.v > 0 or state_->getTargetSpeed() - sd.v < -0.1) {
      pol = MATCH_SPEED;
      delete state_;
      state_ = new MatchSpeed(*this, sd.v, &sd);
      prev_pol_ = pol;
      return;
    // If we just now changed to KEEP_LANE, change the state
    } else if (pol == KEEP_LANE and prev_pol_ != pol) {
      delete state_;
      state_ = new KeepLane(*this, sd.v);
    }
  }

  prev_pol_ = pol;

  // Update current state
  state_->update();
  return;
}
