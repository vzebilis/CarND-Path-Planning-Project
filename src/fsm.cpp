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
// it would take to match target speed with 0 final acceleration
void FSM::matchCarSpeed(double init_s, double init_d, double final_s, double final_d, double final_v)
{
  Movement m = getLastPlannedMovement();
  double v_tot = sqrt(pow(m.vx, 2) + pow(m.vy, 2));
  double a_tot = sqrt(pow(m.ax, 2) + pow(m.ay, 2));
  assert(init_s  >= 0);
  assert(final_v <= MAX_SPEED);

  if (DEBUG) std::cout << "Try to match target speed: " << final_v << std::endl;

  double delta_v = final_v - v_tot;
  // If we are already at the target speed, with zero acceleration, then nothing to do
  if (delta_v > 0 and delta_v < 0.1 and isZero(a_tot)) return;


  // How many timesteps do we need at least to shed our initial acceleration?
  int shed_acc_steps      = ceil((fabs(a_tot) / MAX_JERK) / TIME_RES);
  double delta_speed_shed = 0;
  double a_shed           = 0;
  double a_jerk           = (a_tot > 0)? -MAX_JERK: MAX_JERK;
  // How much speed are we going to add while shedding initial acceleration?
  for (int i = 0; i < shed_acc_steps; ++i) {
    a_shed               += a_jerk * TIME_RES;
    delta_speed_shed     += a_shed * TIME_RES;
  }
  const bool accel = delta_v + delta_speed_shed >= 0;

  // The flip point: speed at which to switch from accelerating to decelerating or vice versa
  double flip_point = v_tot + (delta_v + delta_speed_shed) / 2;
  double T = 0;
  double S = init_s;
  double cur_speed = v_tot;
  double cur_acc = a_tot;
  double cur_x = m.x;
  double cur_yaw = m.yaw;
  while (((accel)? final_v > cur_speed: cur_speed > final_v)) {
    Movement prev_m = m;
    if (DEBUG) {
      std::cout << "Cur speed: " << cur_speed << " Cur acc: " << cur_acc << " Cur S: " << S <<
                   " Speed flip point: " << flip_point << std::endl;
    }
    T           += TIME_RES;
    S           += cur_speed * TIME_RES;

    cur_x += cur_speed * cos(cur_yaw);

    cur_speed   += cur_acc * TIME_RES;
    bool flipped = (accel)? cur_speed > flip_point: cur_speed < flip_point;
    double jerk  = (accel)? MAX_JERK: -MAX_JERK;
    jerk         = (flipped)? -jerk: jerk;
    cur_acc     += jerk * TIME_RES;
    if (fabs(cur_acc) >= MAX_ACC) {
      cur_acc = (cur_acc > 0)? MAX_ACC: -MAX_ACC;
      // Update flip point as well
      flip_point = final_v + ((cur_acc > 0)? -1: 1) * delta_speed_to_zero_acc_;
    }
    if (flipped and isZero(cur_acc)) break;
  }


}

void State::matchTargetSpeed(double d, double trg_v)
{
  if (DEBUG) std::cout << "Target speed to match: " << trg_v << std::endl;
  forward_traj_.clear();
  cur_traj_idx_ = 0;

  Movement prev_m = fsm_.getLastPlannedMovement();
//
//  // Create a local spline for use with local car coords
//  double ref_x, ref_y, ref_yaw;
//  std::unique_ptr<tk::spline> spl(fsm_.createLocalSpline(d, prev_m, ref_x, ref_y, ref_yaw));
//
//  double trg_x = 30; //trg_xy[0];
//  double trg_y = (*spl)(trg_x);
//  double trg_dist = sqrt(pow(trg_x, 2) + pow(trg_y, 2));
//  double N = trg_dist / (TIME_RES * trg_v);
//  double x_add_on = 0; // Always in car-local coords
//  double prev_x_point = prev_m.x; //0;
//  double prev_y_point = prev_m.y; //(*spl)(prev_x_point);
//  double cur_yaw = ref_yaw;


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


  // How many timesteps do we need at least to shed our initial acceleration?
  int shed_acc_steps      = ceil((fabs(tot_a) / MAX_JERK) / TIME_RES);
  double delta_speed_shed = 0;
  double a_shed           = 0;
  double a_jerk           = (tot_a > 0)? -MAX_JERK: MAX_JERK;
  // How much speed are we going to add while shedding initial acceleration?
  for (int i = 0; i < shed_acc_steps; ++i) {
    a_shed               += a_jerk * TIME_RES;
    delta_speed_shed     += a_shed * TIME_RES;
  }
  const bool accel = delta_v + delta_speed_shed >= 0;

  if (DEBUG) std::cout << "MatchTargetSpeed: initial acc: " << tot_a << ". Need to shed for dv: " <<
      delta_speed_shed << std::endl;

  // The flip point: speed at which to switch from accelerating to decelerating or vice versa
  double flip_point = tot_v + (delta_v + delta_speed_shed) / 2;
  if (DEBUG) std::cout << "MatchTargetSpeed: initial flip_point " << flip_point << std::endl;

  //double T = 0;
  double cur_x = prev_m.x;
  double cur_y = prev_m.y;
  double cur_speed = tot_v;
  double cur_acc = tot_a;
  double cur_yaw = prev_m.yaw;
  double prev_x = cur_x;
  double prev_y = cur_y;
  if (DEBUG) std::cout << "MatchTargetSpeed: " << (accel? "accelerating": "decelerating") << std::endl;

  while (((accel)? trg_v > cur_speed: cur_speed > trg_v)) {
    //T           += TIME_RES;
    cur_x       += vx * TIME_RES;
    cur_y       += vy * TIME_RES;
    vx          += ax * TIME_RES;
    vy          += ay * TIME_RES;
    cur_speed   = sqrt(pow(vx, 2) + pow(vy, 2));
    bool flipped = (accel)? cur_speed > flip_point: cur_speed < flip_point;
    if (flipped and DEBUG) std::cout << "FLIPPED!" << std::end;
    double jerk  = (accel)? MAX_JERK: -MAX_JERK;
    jerk         = (flipped)? -jerk: jerk;
    double jx    = jerk * cos(cur_yaw);
    double jy    = jerk * sin(cur_yaw);
    double n_ax  = ax + jx * TIME_RES;
    double n_ay  = ay + jy * TIME_RES;
    cur_acc      = sqrt(pow(n_ax, 2) + pow(n_ay, 2));
    if (fabs(cur_acc) >= MAX_ACC) {
      cur_acc = (cur_acc > 0)? MAX_ACC: -MAX_ACC;
      ax = ((ax > 0)? 1: -1) * cur_acc * cos(cur_yaw);
      ay = ((ay > 0)? 1: -1) * cur_acc * sin(cur_yaw);
      // Update flip point as well
      flip_point = trg_v + (accel? -1: 1) * fsm_.delta_speed_to_zero_acc_;
      if (DEBUG) std::cout << "MatchTargetSpeed: updated flip_point " << flip_point << std::endl;
    } else {
      ax = n_ax;
      ay = n_ay;
    }
    if (!accel) cur_acc = -cur_acc;
    cur_yaw = atan2(cur_y - prev_y, cur_x - prev_x);
    prev_x = cur_x;
    prev_y = cur_y;

    forward_traj_.emplace_back();
    TrajNode & tn = forward_traj_.back();
    tn.x = cur_x;
    tn.y = cur_y;
    auto cur_sd = getFrenet(cur_x, cur_y, cur_yaw, fsm_.map_.x, fsm_.map_.y);
    tn.s = cur_sd[0];
    tn.d = d;

    if (DEBUG) {
      std::cout << "Cur speed: " << cur_speed << " Cur acc: " << cur_acc << " Cur x: " << cur_x <<
          " y: " << cur_y << " vx: " << vx << " vy: " << vy << " ax: " << ax << " ay: " << ay << std::endl;
    }

    if (flipped and isZero(cur_acc)) break;
  }
  if (DEBUG) std::cout << "MatchTargetSpeed: computed trajectory of " << forward_traj_.size() << " steps\n";

//  while (true) {
//    double x_point = x_add_on + trg_x / N;
//    double y_point = (*spl)(x_point);
//    auto world_xy = translateToWorldCoords(x_point, y_point, ref_yaw, ref_x, ref_y);
//    x_add_on = x_point;
//    x_point = world_xy[0];
//    y_point = world_xy[1];
//
//    Movement m = fsm_.getLastPlannedMovement();
//
//    //    auto m = getLastPlannedMovement();
//    //    if (DEBUG) { std::cout << "Prev "; printMovement(prev_m); }
//    if (DEBUG) { std::cout << "Attempt "; printMovement(m); }
//    //    if (updateLimitMovement(m, prev_m, trg_v)) {
//    //      if (DEBUG) { std::cout << "Updated "; printMovement(m); }
//    //      //next_x_vals_.back() = x_point = m.x;
//    //      //next_y_vals_.back() = y_point = m.y;
//    //    }
//
//    // Update cur_yaw
//    cur_yaw = atan2(y_point - prev_y_point, x_point - prev_x_point);
//    prev_x_point = x_point;
//    prev_y_point = y_point;
//
//    // Populate the rest of the values
//    forward_traj_.emplace_back();
//    TrajNode & tn = forward_traj_.back();
//    tn.d = d;
//
//    next_x_vals_.push_back(x_point);
//    next_y_vals_.push_back(y_point);
//    auto cur_sd = getFrenet(x_point, y_point, cur_yaw, map_.x, map_.y);
//    next_s_vals_.push_back(cur_sd[0]);
//    next_d_vals_.push_back(d);
//    prev_acc_.push_back(sqrt(pow(m.ax, 2) + pow(m.ay, 2)));
//    prev_speed_.push_back(sqrt(pow(m.vx, 2) + pow(m.vy, 2)));
//    if (DEBUG) {
//      std::cout << "Added new trajectory step " << (t+1) << " with x: " << x_point << " y: " << y_point <<
//          " keep_speed: 1" << std::endl;
//    }
//    prev_m = m;
//  }
}

// Method to compute and return the minimum amount of road length and time
// it would take to match target speed with 0 final acceleration
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
    if (fabs(cur_acc) >= MAX_ACC) {
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
  ret.time        = T;

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
void State::computeXYTrajectory(TrajData & td)
{
  assert(td.time);
  Context & cxt = *fsm_.cxt_;
  size_t sz = cxt.path_x.size();
  if (sz +  cur_traj_idx_ > forward_traj_.size()) forward_traj_.clear();
  else forward_traj_.erase(forward_traj_.begin(), forward_traj_.begin() + sz + cur_traj_idx_);

  size_t time_steps = ceil(td.time / TIME_RES);

  // Here we need to start with two simple points close to our reference location
  vector<double> xs, ys;
  double ref_yaw = deg2rad(cxt.yaw);
  double ref_x   = cxt.x;
  double ref_y   = cxt.y;
  double ref_v   = td.init_speed;
  double trg_v   = td.final_speed;
  bool accel     = trg_v > ref_v;
  if (sz < 2) {
    double prev_x = ref_x - cos(ref_yaw);
    double prev_y = ref_y - sin(ref_yaw);
    xs.push_back(prev_x); ys.push_back(prev_y);
  } else {
    ref_x = cxt.path_x[sz-1];
    ref_y = cxt.path_y[sz-1];
    double ref_x_prev = cxt.path_x[sz-2];
    double ref_y_prev = cxt.path_y[sz-2];
    ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);
    xs.push_back(ref_x_prev); ys.push_back(ref_y_prev);
  }
  xs.push_back(ref_x); ys.push_back(ref_y);

  // Now we will add some more points further away
  auto next_wp0 = getXY(cxt.s + 30, td.final_d, fsm_.map_.s, fsm_.map_.x, fsm_.map_.y);
  auto next_wp1 = getXY(cxt.s + 60, td.final_d, fsm_.map_.s, fsm_.map_.x, fsm_.map_.y);
  auto next_wp2 = getXY(cxt.s + 90, td.final_d, fsm_.map_.s, fsm_.map_.x, fsm_.map_.y);
  xs.push_back(next_wp0[0]); ys.push_back(next_wp0[1]);
  xs.push_back(next_wp1[0]); ys.push_back(next_wp1[1]);
  xs.push_back(next_wp2[0]); ys.push_back(next_wp2[1]);

  // Now we need to transform these points to the car-local coordinates:
  // The reference point is at (0,0) and the x axis points to our yaw angle.
  // This way we can use the spline as f(x), and it would still behave as a function (single value of y for x)
  for (int i = 0; i < xs.size(); ++i) {
    double shift_x = xs[i] - ref_x;
    double shift_y = ys[i] - ref_y;
    xs[i] = (shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw));
    ys[i] = (shift_x * sin(-ref_yaw) + shift_y * cos(-ref_yaw));
  }

  tk::spline spl;
  spl.set_points(xs, ys);

  // How many point do I need to sample from the spline so that we still keep under our
  // velocity and acceleration constraints?
  auto trg_xy = getXY(td.final_s, td.final_d, fsm_.map_.s, fsm_.map_.x, fsm_.map_.y);
  double trg_x = trg_xy[0];
  double trg_y = spl(trg_x);
  double trg_dist = sqrt(pow(trg_x, 2) + pow(trg_y, 2));
  double N = time_steps; //trg_dist / (TIME_RES * ref_v);
  double x_add_on = 0;
  double prev_x_point = 0;
  double prev_y_point = spl(prev_x_point);
  double cur_yaw = ref_yaw;
  for (size_t t = 0; t < time_steps; ++t) {
    double x_point = x_add_on + trg_x / N;
    double y_point = spl(x_point);
    x_add_on = x_point;
    double x_ref = x_point;
    double y_ref = y_point;
    // Re-translate to original world coordinates
    x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
    y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

    x_point += ref_x;
    y_point += ref_y;

    // Update cur_yaw
    cur_yaw = atan2(y_point - prev_y_point, x_point - prev_x_point);
    prev_x_point = x_point;
    prev_y_point = y_point;

    TrajNode tn; tn.x = x_point; tn.y = y_point; tn.s = 0; tn.d = 0;
    forward_traj_.push_back(tn);
//    fsm_.next_x_vals_.push_back(x_point);
//    fsm_.next_y_vals_.push_back(y_point);
  }

  return;
}
//
// Generic state method to compute the necessary trajectory based on
// the supplied TrajData info
//
void State::computeTrajectory(TrajData & td)
{
  assert(td.time);
//  size_t sz = cur_traj_idx_ + fsm_.cxt_->path_x.size();
//  if (sz > forward_traj_.size()) forward_traj_.clear();
//  else forward_traj_.erase(forward_traj_.begin() + sz, forward_traj_.end());
  forward_traj_.clear();
  cur_traj_idx_  = 0;


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

//
// Function to check which of the available lanes is best to use
//
SensorData FSM::checkLanes(const PositionData & p)
{
  if (DEBUG) std::cout << "checkLanes: current lane: " << printLane(lane_) << std::endl;
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
  // DEBUG DEBUG DEBUG DEBUG
  // DEBUG DEBUG DEBUG DEBUG
  // DEBUG DEBUG DEBUG DEBUG
  // TODO: for now skip changing lanes
  return best;

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
  //TrajData td = computeTargetPos(p_, sd);
  fsm_.prev_trg_speed_ = trg_speed_;
  // TODO: hack!
  p_.d = sd->d;
  //computeTrajectory(td); // Base class method
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
  //p_ = fsm_.getLastProcessed();
  //fsm_.clearAll();
  p_ = fsm_.getMostRecentNotProcessed();
  fsm_.clearProcessed();
  //TrajData td = computeTargetPos(p_, sd);
  matchTargetSpeed(sd? sd->d: p_.d, trg_speed);
  fsm_.prev_trg_speed_ = trg_speed_;
  //computeXYTrajectory(td);
  //computeTrajectory(td); // Base class method
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
  //fsm_.incrementByProcessed(cur_traj_idx_);
  // Apply the new trajectory by creating the new next_x_vals and next_y_vals
  fsm_.applyTrajectory(p_, forward_traj_, cur_traj_idx_);
}

static TrajData destTrajData(PositionData & p, SensorData * sd)
{
  TrajData td;
  td.time = 3.5;

  td.init_s = p.s;
  td.final_s = p.s + 30;

  td.init_d = p.d;
  td.final_d = p.d;

  td.init_speed = p.v;
  td.final_speed = MAX_SPEED * MAX_SPEED_MAR;

  td.init_acc = p.a;
  td.final_acc = 0;

  if (!sd or sd->s == -1) {
    if (DEBUG) std::cout << "Setting target speed to MAX\n";
    return td;
  }

  double trg_speed = sd->v;
  trg_speed        = fmin(trg_speed, MAX_SPEED * MAX_SPEED_MAR);
  td.final_speed   = trg_speed;

  // Don't waste time
  if (trg_speed >= p.v) {
    if (DEBUG) std::cout << "Setting target speed to " << trg_speed << std::endl;
    return td;
  }

  // Otherwise, decelerate more smoothly if we have enough space/time to do it
  // In order to so so, we calculate the future position of the observed car
  // in the time it will take us to match the target speed.
  // This give us an extra margin to decelerate more smoothly
  double car_s = sd->s;
  if (car_s < td.init_s) car_s += MAX_S;
  if (car_s - FRONT_CAR_MAR > td.init_s) car_s -= FRONT_CAR_MAR;
  auto time_trav = getTimeToIntercept(car_s - td.init_s, td.init_speed, trg_speed);
  if (DEBUG) std::cout << "Time to intercept car: " << time_trav.first << " and s: " << time_trav.second << std::endl;
  td.final_s     = td.init_s + time_trav.second;
  td.final_speed = trg_speed;
  td.time        = time_trav.first;

  return td;
}

//
// Compute target position, for matching observed speed
//
TrajData MatchSpeed::computeTargetPos(PositionData & p, SensorData * sd)
{
  //return destTrajData(p, sd);
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
  td.init_d = td.final_d = init_d;
  return td;
}

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
    //prev_x = cxt_->x - cos(cxt_->yaw);
    //prev_y = cxt_->y - sin(cxt_->yaw);
  } else {
    double x1 = next_x_vals_[sz-1];
    double y1 = next_y_vals_[sz-1];
    double x2 = next_x_vals_[sz-2];
    double y2 = next_y_vals_[sz-2];
    double x3 = next_x_vals_[sz-3];
    double y3 = next_y_vals_[sz-3];
    // TODO: assuming that is x1 == x2 then we are a right angle to the x axis.
    if (x1 == x2) {
      m.yaw = cxt_->yaw;//((y1 - y2 > 0)? 1 : -1 ) * pi() / 2;
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
 * Helper function to limit the proposed movement so that no constraint violations happen
 * @param m Proposed movement
 * @param prev_m Previous movement
 * @param trg_v Target velocity to use as a constraint
 * @return True if there was an limiting update, False otherwise
 */
static bool updateLimitMovement(Movement & m, const Movement & prev_m, double trg_v)
{
  double & cur_vx = m.vx; //prev_speed * cos(yaw);
  double & cur_vy = m.vy; //prev_speed * sin(yaw);
  double & cur_ax = m.ax; //cur_vx - prev_vx;
  double & cur_ay = m.ay; //cur_vy - prev_vy;
  double prev_vx = prev_m.vx;
  double prev_vy = prev_m.vy;
  double prev_ax = prev_m.ax;
  double prev_ay = prev_m.ay;
//  double prev_tot_v = sqrt(pow(prev_vx, 2) + pow(prev_vy, 2));
//  double prev_tot_a = sqrt(pow(prev_ax, 2) + pow(prev_ay, 2));
//  double lim_vx = MAX_SPEED * MAX_SPEED_MAR * cos(m.yaw);
//  double lim_vy = MAX_SPEED * MAX_SPEED_MAR * sin(m.yaw);
//  double lim_ax = MAX_ACC * cos(m.yaw);
//  double lim_ay = MAX_ACC * sin(m.yaw);
//  double lim_jx = MAX_JERK * cos(m.yaw);
//  double lim_jy = MAX_JERK * sin(m.yaw);

  double lim_v = trg_v;
  double lim_a = MAX_ACC;
  double lim_j = MAX_JERK;

  bool updated = false;

  //
  // Limit velocity to the target one
  //
  double tot_v = sqrt(pow(cur_vx, 2) + pow(cur_vy, 2));
  if (tot_v > lim_v) {
    updated = true;
    if (DEBUG) std::cout << "limiting speed " << tot_v;
    tot_v = lim_v;
    cur_vx = tot_v * cos(prev_m.yaw) * (cur_vx > 0? 1: -1);
    cur_vy = tot_v * sin(prev_m.yaw) * (cur_vy > 0? 1: -1);
    if (DEBUG) std::cout << " tot_v: " << tot_v << " cur_vx: " << cur_vx << " cur_vy: " << cur_vy << std::endl;
    cur_ax = (cur_vx - prev_vx) / TIME_RES;
    cur_ay = (cur_vy - prev_vy) / TIME_RES;
  }
  //double tot_a = tot_v - prev_tot_v;

  //
  // Limit acceleration to +/- max one and update velocity accordingly
  //
  double tot_a = sqrt(pow(cur_ax, 2) + pow(cur_ay, 2));
  if (tot_a > lim_a or tot_a < -lim_a) {
    updated = true;
    if (DEBUG) std::cout << "limiting acceleration " << tot_a;
    tot_a = (tot_a > 0)? lim_a: -lim_a;
    cur_ax = tot_a * cos(prev_m.yaw) * (cur_ax > 0? 1: -1);
    cur_ay = tot_a * sin(prev_m.yaw) * (cur_ay > 0? 1: -1);
    if (DEBUG) std::cout << " tot_a: " << tot_a << " cur_ax: " << cur_ax << " cur_ay: " << cur_ay << std::endl;
    cur_vx = prev_vx + cur_ax * TIME_RES;
    cur_vy = prev_vy + cur_ay * TIME_RES;
  }

  //
  // Limit jerk to +/- max one and update acceleration and velocity accordingly
  //
  double cur_jx = (cur_ax - prev_ax) / TIME_RES;
  double cur_jy = (cur_ay - prev_ay) / TIME_RES;
  //double tot_j = tot_a - prev_tot_a;
  double tot_j = sqrt(pow(cur_jx, 2) + pow(cur_jy, 2));
  if (tot_j > lim_j or tot_j < -lim_j) {
    updated = true;
    if (DEBUG) std::cout << "limiting jerk " << tot_j;
    tot_j = (tot_j > 0)? lim_j: -lim_j;
    cur_jx = tot_j * cos(prev_m.yaw) * (cur_jx > 0? 1: -1);
    cur_jy = tot_j * sin(prev_m.yaw) * (cur_jy > 0? 1: -1);
    if (DEBUG) std::cout << " tot_j: " << tot_j << " cur_jx: " << cur_jx << " cur_jy: " << cur_jy << std::endl;
    cur_ax = prev_ax + cur_jx * TIME_RES;
    cur_ay = prev_ay + cur_jy * TIME_RES;
    cur_vx = prev_vx + cur_ax * TIME_RES;
    cur_vy = prev_vy + cur_ay * TIME_RES;
  }
  // Update the movement X, Y coords if VX and VY actually changed
  if (updated) {
    m.x = prev_m.x + m.vx * TIME_RES;
    m.y = prev_m.y + m.vy * TIME_RES;
  }

  return updated;
//
//  if (cur_vx > lim_vx or cur_vx < -lim_vx) {
//    cur_vx = (cur_vx > 0)? lim_vx: -lim_vx;
//    cur_ax = (cur_vx - prev_vx) / TIME_RES;
//  }
//  if (cur_vy > lim_vy or cur_vy < -lim_vy) {
//    cur_vy = (cur_vy > 0)? lim_vy: -lim_vy;
//    cur_ay = (cur_vy - prev_vy) / TIME_RES;
//  }
//  if (cur_ax > lim_ax or cur_ax < -lim_ax) {
//    cur_ax = (cur_ax > 0)? lim_ax: -lim_ax;
//    cur_vx = prev_vx + cur_ax * TIME_RES;
//  }
//  if (cur_ay > lim_ay or cur_ay < -lim_ay) {
//    cur_ay = (cur_ay > 0)? lim_ay: -lim_ay;
//    cur_vy = prev_vy + cur_ay * TIME_RES;
//  }
//
//  double cur_jx = (cur_ax - prev_ax) / TIME_RES;
//  if (cur_jx > MAX_JERK or cur_jx < -MAX_JERK) {
//    cur_jx = (cur_jx > 0)? MAX_JERK: -MAX_JERK;
//    cur_ax = prev_ax + cur_jx * TIME_RES;
//    cur_vx = prev_vx + cur_ax * TIME_RES;
//  }
//  double cur_jy = (cur_ay - prev_ay) / TIME_RES;
//  if (cur_jy > MAX_JERK or cur_jy < -MAX_JERK) {
//    cur_jy = (cur_jy > 0)? MAX_JERK: -MAX_JERK;
//    cur_ay = prev_ay + cur_jy * TIME_RES;
//    cur_vy = prev_vy + cur_ay * TIME_RES;
//  }
}

/**
 * A method to compute the next world X and Y coordinates and add them to the FSM.
 * @param start_x The starting point in the x axis in car-local coords.
 * @param step_x The step to try and use. Will be limited if violating constraints.
 * @param spl The Spline to use for interpolation
 * @param ref_yaw The reference point world YAW for translation between car and world coords.
 * @param ref_x The reference point world X for translation between car and world coords.
 * @param ref_y The reference point world Y for translation between car and world coords.
 */
void FSM::computeNextXY(double start_x, double step_x, tk::spline & spl,
    double ref_yaw, double ref_x, double ref_y)
{
  double x_point = start_x + step_x;
  double y_point = spl(x_point);

  double x_ref = x_point;
  double y_ref = y_point;

  // Re-translate to original world coordinates
  x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
  y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));
  x_point += ref_x;
  y_point += ref_y;

  next_x_vals_.push_back(x_point);
  next_y_vals_.push_back(y_point);
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
//  size_t sz = cxt.path_x.size();
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
//    ref_x = cxt.path_x[sz-1];
//    ref_y = cxt.path_y[sz-1];
//    double ref_x_prev = cxt.path_x[sz-2];
//    double ref_y_prev = cxt.path_y[sz-2];
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

  double d0, d1, d2;
  d0 = d1 = d2 = d;

//  if (fabs(cur_d - d) > 3) {
//    d0 += (cur_d - d0) / 2;
//  }

  // Now we will add some more points further away
  auto next_wp0 = getXY(cur_sd[0] + 30, d0, map_.s, map_.x, map_.y);
  auto next_wp1 = getXY(cur_sd[0] + 60, d1, map_.s, map_.x, map_.y);
  auto next_wp2 = getXY(cur_sd[0] + 90, d2, map_.s, map_.x, map_.y);
  xs.push_back(next_wp0[0]); ys.push_back(next_wp0[1]);
  xs.push_back(next_wp1[0]); ys.push_back(next_wp1[1]);
  xs.push_back(next_wp2[0]); ys.push_back(next_wp2[1]);

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
 * @param trg_v Traget speed
 * @param d Current d Frenet value
 */
void FSM::createXYFromSpline(double trg_v, double d,  Movement & prev_m)
{
  // TODO: for now overriding to the previous trg_speed of the FSM
  //trg_v = prev_trg_speed_;
  if (DEBUG) std::cout << "Target speed for Spline: " << trg_v << std::endl;

  size_t sz = next_x_vals_.size();

  // Create a local spline for use with local car coords
  double ref_x, ref_y, ref_yaw;
  std::unique_ptr<tk::spline> spl(createLocalSpline(d, prev_m, ref_x, ref_y, ref_yaw));

  // How many point do I need to sample from the spline so that we still keep under our
  // velocity and acceleration constraints?
//  auto trg_xy = getXY(td.final_s, td.final_d, map_.s, map_.x, map_.y);
  double trg_x = 30; //trg_xy[0];
  double trg_y = (*spl)(trg_x);
  double trg_dist = sqrt(pow(trg_x, 2) + pow(trg_y, 2));
  double N = trg_dist / (TIME_RES * trg_v);
  double x_add_on = 0; // Always in car-local coords
  double prev_x_point = prev_m.x; //0;
  double prev_y_point = prev_m.y; //(*spl)(prev_x_point);
  double cur_yaw = ref_yaw;
  for (size_t t = sz; t < TRAJ_STEPS; ++t) {
    double x_point = x_add_on + trg_x / N;
    double y_point = (*spl)(x_point);
    auto world_xy = translateToWorldCoords(x_point, y_point, ref_yaw, ref_x, ref_y);
    x_add_on = x_point;
    x_point = world_xy[0];
    y_point = world_xy[1];

    Movement m = getLastPlannedMovement();

//    auto m = getLastPlannedMovement();
//    if (DEBUG) { std::cout << "Prev "; printMovement(prev_m); }
    if (DEBUG) { std::cout << "Attempt "; printMovement(m); }
//    if (updateLimitMovement(m, prev_m, trg_v)) {
//      if (DEBUG) { std::cout << "Updated "; printMovement(m); }
//      //next_x_vals_.back() = x_point = m.x;
//      //next_y_vals_.back() = y_point = m.y;
//    }

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
    prev_acc_.push_back(sqrt(pow(m.ax, 2) + pow(m.ay, 2)));
    prev_speed_.push_back(sqrt(pow(m.vx, 2) + pow(m.vy, 2)));
    if (DEBUG) {
      std::cout << "Added new trajectory step " << (t+1) << " with x: " << x_point << " y: " << y_point <<
          " keep_speed: 1" << std::endl;
    }
    prev_m = m;
  }
}


//
// Method to apply the pre-computed trajectory
//
void FSM::applyTrajectory(PositionData & p, vector<TrajNode> & forward_traj, size_t & cur_traj_idx)
{
  double prev_ns   = p.s;
  double prev_d    = p.d;
  double cur_speed = p.v;
  double cur_acc   = p.a;

  double prev_speed = cur_speed;

  updateLane(prev_d);
  static Movement prev_m;
  prev_m = getLastPlannedMovement();

  bool keep_speed = (cur_traj_idx >= forward_traj.size());

  size_t i = next_x_vals_.size();

  if (keep_speed) {
//    auto m = getLastPlannedMovement();
    createXYFromSpline(prev_trg_speed_, prev_d, prev_m);
    return;
  }

  if (DEBUG) {
    printPositionData("ApplyTrajectory", p);
    std::cout << "cur_traj_idx: " << cur_traj_idx << " out of " << forward_traj.size() << std::endl;
  }

  double time = 0;
  // Build the rest of the steps

  for (; i < TRAJ_STEPS; ++i) {
    // If we completed our computed trajectory, just change policy to keep current speed
    if (!keep_speed and cur_traj_idx >= forward_traj.size()) {
      //TrajData td = initTrajData(i, prev_ns, prev_speed, prev_d);
      //spd = getShortSplines(td);
      keep_speed = true;
      time = 0;
      if (DEBUG) std::cout << "Finished generated trajectory. Now keeping speed (TrajSize: " << forward_traj.size() << ")\n";
      createXYFromSpline(prev_trg_speed_, prev_d, prev_m);
      return;
    }
    time += TIME_RES;
    double nx, ny, ns, nd;
    if (keep_speed) {
//      if (DEBUG) std::cout << "PrevYaw: " << prev_m.yaw << " cos: " << cos(prev_m.yaw) << " sin: " << sin(prev_m.yaw) <<std::endl;
//      //double prev_vx = prev_speed * cos(yaw);
//      //double prev_vy = prev_speed * sin(yaw);
//      //double s_step = sqrt(pow(prev_vx * TIME_RES, 2) + pow(prev_vy * TIME_RES, 2));
//      auto prev_xy = getXY(prev_ns, prev_d, map_.s, map_.x, map_.y);
//
//      //nx = prev_xy[0] + m.vx * TIME_RES;
//      //ny = prev_xy[1] + m.vy * TIME_RES;
//      nx = next_x_vals_.back() + prev_m.vx * TIME_RES;
//      ny = next_y_vals_.back() + prev_m.vy * TIME_RES;
//      //double s_step = fmax(prev_speed * cos(yaw), prev_speed * sin(yaw)) * TIME_RES;
//      //ns = prev_ns + s_step;//prev_speed * TIME_RES;
//      auto next_sd = getFrenet(nx, ny, prev_m.yaw, map_.x, map_.y);
//      ns = next_sd[0];
//      nd = prev_d;
//      //ns = fmod(ns, MAX_S);
////      nd = prev_d;
////      auto vecXY = getXYfromSpline(ns, nd);
////      nx = vecXY[0];
////      ny = vecXY[1];
    } else {
      nx = forward_traj[cur_traj_idx].x;
      ny = forward_traj[cur_traj_idx].y;
      ns = forward_traj[cur_traj_idx].s;
      nd = forward_traj[cur_traj_idx].d;
    }
    ++cur_traj_idx;
    next_x_vals_.push_back(nx);
    next_y_vals_.push_back(ny);
//    auto m = getLastPlannedMovement();
//    if (DEBUG) { std::cout << "Prev "; printMovement(prev_m); }
//    if (DEBUG) { std::cout << "Attempted "; printMovement(m); }
//    if (updateLimitMovement(m, prev_m, prev_trg_speed_)) {
//      if (DEBUG) { std::cout << "Updated "; printMovement(m); }
//      next_x_vals_.back() = nx = m.x;
//      next_y_vals_.back() = ny = m.y;
//      double cur_yaw = atan2(m.y - prev_m.y, m.x - prev_m.x);
//      auto cur_sd = getFrenet(m.x, m.y, cur_yaw, map_.x, map_.y);
//      ns = cur_sd[0];
//    }
    next_s_vals_.push_back(ns);
    next_d_vals_.push_back(nd);

    prev_m     = getLastPlannedMovement();//m;

    if (DEBUG) {
      std::cout << "Added new trajectory step " << (i+1) << " with x: " << nx << " y: " << ny << " s: " <<
          ns << " d: " << nd << " speed: " << cur_speed << " acc: " << cur_acc << " keep_speed: " <<
          keep_speed << std::endl << std::endl;
    }

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
  prev_pol_ = MATCH_SPEED;
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
    std::cout << "[Step " << stepCnt << "] Current internal state  (s,d,x,y,yaw,v,a): (" <<
        p.s << ", " << p.d << ", " << p.x << ", " << p.y << ", " << deg2rad(cxt_->yaw) << ", " << p.v <<
        ", " << p.a << ")\n";
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
    } else if (state_->getTargetSpeed() - sd.v > 0.1 or state_->getTargetSpeed() - sd.v < -0.1) {
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
