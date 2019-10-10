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

  double prev_x = spd_.sx(map_.s[prev_wp]);
  double prev_y = spd_.sy(map_.s[prev_wp]);
  double x2 = spd_.sx(map_.s[wp2]);
  double y2 = spd_.sy(map_.s[wp2]);
  double heading = atan2(y2 - prev_y, x2 - prev_x);

  double seg_x  = spd_.sx(s);
  double seg_y  = spd_.sy(s);
  double seg_dx = spd_.sdx(s);
  double seg_dy = spd_.sdy(s);
  double d = sqrt(pow(seg_dx, 2) + pow(seg_dy, 2)) * d_offset; // Middle lane

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};
}

FSM::FSM(const Map & map) : 
  policy_(MATCH_SPEED), prev_acc_(TRAJ_STEPS, 0), prev_speed_(TRAJ_STEPS, 0), map_(map), prev_trg_speed_(0),
  cur_traj_idx_(0), delta_speed_to_zero_acc_(0)
{
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
  spd_.sx.set_points(map_.s, map_.x);
  spd_.sy.set_points(map_.s, map_.y);
  spd_.sdx.set_points(map_.s, map_.dx);
  spd_.sdy.set_points(map_.s, map_.dy);
}

// Method to compute and return the minimum amount of road length and time
// it would take to acceleratate to target speed with 0 final acceleration
TrajData FSM::computeMatchTargetSpeed(double init_s, double init_speed, double init_acc, double trg_speed) 
{
  assert(init_s     >= 0);
  assert(init_speed >= 0);
  assert(trg_speed  <= MAX_SPEED);

  if (DEBUG) std::cout << "Try to match target speed: " << trg_speed << std::endl;

  TrajData ret; 
  ret.final_speed         = init_speed; // default return value

  //double delta_speed      = trg_speed - init_speed;
  //bool accel              = delta_speed > 0;
  //delta_speed             = fabs(delta_speed);

  // If we are already at the target speed, with zero acceleration, then nothing to do
  if (isZero(trg_speed - init_speed) and isZero(init_acc)) return ret;


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
  // Update the desired delta_speed
  //delta_speed            += delta_speed_shed;
  //if (delta_speed < 0) accel = false;
  const bool accel = trg_speed - init_speed + delta_speed_shed >= 0;

  // The flip point: speed at which to switch from accelerating to decelerating or vice versa
  //double flip_point = delta_speed / 2;
  double flip_point = init_speed + (trg_speed - init_speed + delta_speed_shed) / 2;
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
    //if (cur_speed > MAX_SPEED) cur_speed = MAX_SPEED;
    //else if (cur_speed < 0)    cur_speed = 0;
    //double dv    = cur_acc * TIME_RES;
    //if (cur_speed + dv > MAX_SPEED) dv = MAX_SPEED - cur_speed;
    //else if (cur_speed + dv < 0)    dv = -cur_speed;
    //delta_speed -= dv;
    //cur_speed   += dv;
    //jerk         = (delta_speed > flip_point)? jerk: -jerk;
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
  ret.time        = T;

  return ret;
}

// Method to compute and return the minimum amount of road length and time
// it would take to acceleratate to target speed with 0 final acceleration
//TrajData FSM::computeAccelerateToTrgSpeed(double init_s, double init_speed, double init_acc, double trg_speed) 
//{
//  // How many timesteps do we need at least to shed our initial acceleration?
//  int shed_acc_steps      = ceil((fabs(init_acc) / MAX_JERK) / TIME_RES);
//  double delta_speed      = trg_speed - init_speed;
//  double delta_speed_shed = 0;
//  double a_shed           = 0;
//  // How much speed are we going to add while shedding initial acceleration?
//  for (int i = 0; i < shed_acc_steps; ++i) {
//    a_shed += MAX_JERK * TIME_RES;
//    delta_speed_shed += a_shed * TIME_RES;
//  }
//  // Update the desired delta_speed
//  delta_speed -= delta_speed_shed;
//  // The flip point, to go from accelerating to decelerating
//  // is the new delta_speed / 2
//  double flip_point = (delta_speed > 0)? delta_speed / 2: 0;
//  double T = 0;
//  double S = init_s;
//  double cur_speed = init_speed;
//  double cur_acc = init_acc;
//  while (cur_speed < trg_speed) {
//    std::cout << "Cur speed: " << cur_speed << " Cur acc: " << cur_acc << " Cur S: " << S <<
//                 " Flip point: " << flip_point << std::endl;
//    T          += TIME_RES;
//    S          += cur_speed * TIME_RES;
//    cur_speed  += cur_acc * TIME_RES;
//    if (cur_speed > MAX_SPEED) cur_speed = MAX_SPEED;
//    double jerk = ((trg_speed - cur_speed) > flip_point)? MAX_JERK: -MAX_JERK;
//    cur_acc    += jerk * TIME_RES;
//    if (cur_acc > MAX_ACC) {
//      cur_acc = MAX_ACC;
//      // Update flip point as well
//      flip_point = delta_speed_to_zero_acc_;
//    }
//    if (cur_acc < 0) {
//      cur_acc = 0;
//      break;
//    }
//  }
//
//  TrajData ret;
//  ret.init_s      = init_s;
//  ret.init_speed  = init_speed;
//  ret.init_acc    = init_acc;
//  ret.final_s     = S;
//  ret.final_speed = trg_speed;
//  ret.final_acc   = cur_acc;
//  ret.time        = T;
//
//  return ret;
//}

void FSM::computeTrajectory(TrajData & traj_data)
{
  forward_traj_.clear();
  cur_traj_idx_ = 0;
  auto coeffs = computePolyCoefficients(traj_data);
  // Generate trajectory until we reach the required time
  double time = TIME_RES;
  while (time < traj_data.time) {
    forward_traj_.emplace_back();
    auto & ft = forward_traj_.back();
    double new_s = getBestStep(coeffs, time);
    new_s = fmod(new_s, MAX_S);
    auto vecXY = getXYfromSpline(new_s, 6);
    ft.x = vecXY[0];
    ft.y = vecXY[1];
    ft.s = new_s;
    time += TIME_RES;
  }
}

static SensorData getSensorFusion(State & st)
{
  // Get data for car closest to us in s, going forward
  int best_idx  = -1;
  double best_s = MAX_S;
  for (int i = 0; i < st.sensor_fusion.size(); ++i) {
    // Only interested about cars in the same lane
    double cur_d = st.sensor_fusion[i][6]; // get the d value
    if (cur_d < st.d - 2 or cur_d > st.d + 2) continue;

    // Only interested about cars in a reasonable distance
    double cur_s = st.sensor_fusion[i][5]; // get the s value
    if (cur_s < st.s) cur_s += MAX_S; // wraparound
    if (cur_s - st.s > 50) continue;

    if (cur_s - st.s < best_s) {
      best_s   = cur_s - st.s;
      best_idx = i;
    }
  }
  SensorData sd; 
  // Are there any cars detectable in front of us?
  if (best_idx == -1 ) {
    sd.s = -1;
    return sd;
  }
  sd.id = st.sensor_fusion[best_idx][0];
  sd.x  = st.sensor_fusion[best_idx][1];
  sd.y  = st.sensor_fusion[best_idx][2];
  sd.vx = st.sensor_fusion[best_idx][3];
  sd.vy = st.sensor_fusion[best_idx][4];
  sd.s  = st.sensor_fusion[best_idx][5];
  sd.d  = st.sensor_fusion[best_idx][6];
  return sd;
}

// Compute target speed, based on two cases
// i.  if we are currently keeping a steady speed
// ii. if we are still in the middle of a maneuver
TrajData FSM::computeTargetSpeed(double init_s, double init_speed, double init_acc, SensorData & sd)
{
  // Default TrajData to return
  TrajData td;

  // If no car spotted, go to max speed
  if (sd.s == -1) {
    if (DEBUG) std::cout << "Setting target speed back to MAX\n";
    return computeMatchTargetSpeed(init_s, init_speed, init_acc, MAX_SPEED * MAX_SPEED_MAR);
  } 

  // Compute target speed
  double car_s     = sd.s;
  if (car_s < init_s) car_s += MAX_S;
  double trg_speed = sqrt(pow(sd.vx, 2) + pow(sd.vy, 2)) * MPH2MPS;
  trg_speed        = fmin(trg_speed, MAX_SPEED * MAX_SPEED_MAR);
  double delta_v   = trg_speed - init_speed;

  // If we are already close to the target speed, with a steady speed, do not alter course
  //if (delta_v > 0 and delta_v < 2) {
  //  if (DEBUG) std::cout << "Not altering initial speed due to it being close to the target one\n";
  //  td.final_speed = init_speed;
  //  return td; 
  //}

  // Compute a trial trajectory to see at what s and time we will have matched car speed
  td = computeMatchTargetSpeed(init_s, init_speed, init_acc, trg_speed);
  // If we don't have enought time to stop, violate acceleration constraints
  double car_future_s = car_s + td.time * trg_speed;
  // Keep a safety margin from the car in front
  car_future_s -= FRONT_CAR_MAR;

  // Is there is time to reach target speed before the front car future location?
  // if so, do the matching more gradually by using the mean speed between the two

  double trg_s;
  // Decelerate
  if (trg_speed < init_speed) {
    trg_s = car_future_s;
    td.time += (car_future_s - td.final_s) / init_speed;
  } 
  // Accelerate
  if (trg_speed > init_speed) {
    assert(td.final_s < car_future_s); // Sanity check
    // If to MAX_SPEED or close enough, just accelerate to it
    if (trg_speed >= MAX_SPEED * MAX_SPEED_MAR) {
      trg_s = td.final_s;
    // Otherwise try to see if accelerating past it is possible
    } else {
      //auto td2 = computeMatchTargetSpeed(init_s, init_speed, init_acc, trg_speed);
      // TODO: for now keeping conservative speed
      trg_s = td.final_s;
    }
  }

  td.final_s = trg_s;
  return td;
}

//
// Main update function for the FSM based on new State
//
void FSM::update(State & st) 
{
  int sz = st.path_x.size();

  // How many points left in the previous path?
  if (sz == TRAJ_STEPS) {
    if (DEBUG) std::cout << "Returning same lists since nothing was processed\n";
    return; // nothing processed, the lists remain the same
  }

  static double trg_speed = MAX_SPEED * MAX_SPEED_MAR; // initialize

  SensorData sd;
  // Get sensor data
  if (sz) {
    sd = getSensorFusion(st);
    if (sd.s != -1) {
      if (DEBUG) {
        double tot_speed  = sqrt(pow(sd.vx, 2) + pow(sd.vy, 2)) * MPH2MPS;
        std::cout << "Spotted closest car ahead in the same lane at S: " << sd.s << " D: " << sd.d << 
          " VX: " << sd.vx << " VY: " << sd.vy << " TotalSpeed: " << tot_speed << std::endl;
      }
    } else {
      if (DEBUG) std::cout << "No car spotted in the same lane!\n";
    }
  }



  int processed = TRAJ_STEPS - sz;

  double last_x_processed     = sz? next_x_vals_[processed -1]: st.x;
  double last_y_processed     = sz? next_y_vals_[processed -1]: st.y;
  double last_s_processed     = sz? next_s_vals_[processed -1]: st.s;
  double last_speed_processed = sz? prev_speed_[processed -1]: 0;
  double last_acc_processed   = sz? prev_acc_[processed -1]: 0;

  // First copy over the old trajectory
  if (sz) {
    next_x_vals_.clear();
    next_y_vals_.clear();
    next_x_vals_.insert(next_x_vals_.begin(), st.path_x.begin(), st.path_x.end());
    next_y_vals_.insert(next_y_vals_.begin(), st.path_y.begin(), st.path_y.end());
    // Erase processed first elements 
    next_s_vals_.erase(next_s_vals_.begin(), next_s_vals_.begin() + processed);
    prev_acc_.erase(prev_acc_.begin(), prev_acc_.begin() + processed);
    prev_speed_.erase(prev_speed_.begin(), prev_speed_.begin() + processed);
  }

  if (DEBUG) {
    std::cout << "=======================================================================\n";
    std::cout << "Processed: " << (sz? std::to_string(processed): "N/A") << std::endl;
  }

  if (DEBUG) {
    static int stepCnt = 0;
    static double last_st_s     = st.s;
    static double last_st_speed = 0;
    double st_end_s   = (sz)? st.end_s: st.s;
    double st_speed   = ((st.s - last_st_s) / TIME_RES) / (sz? processed: 1); 
    double st_speed_i = (sz)? st.speed * MPH2MPS: 0;  // instant speed
    double st_acc     = ((st_speed_i - last_st_speed) / TIME_RES) / (sz? processed: 1);
    last_st_speed     = st_speed_i;
    last_st_s         = st.s;
    stepCnt += (sz)? processed: 0;
    std::cout << "[Step " << stepCnt << "] Current given state (s,end_s,x,y,speed,speed_i,acc): (" << st.s << 
      ", " << st_end_s << ", " << st.x << ", " << st.y << ", " << st_speed << ", " << st_speed_i <<
      ", " << st_acc << ")\n";
    if (sz) {
      std::cout << "[Step " << stepCnt << "] Current internal state          (s,x,y,speed,acc): (" << 
        last_s_processed << ",        " << last_x_processed << ", " << last_y_processed << ", " << 
        last_speed_processed << ", " << last_acc_processed << ")\n";
    }
  }

  TrajData td;
  // See if we need to adjust trajectory
  if (sz) {
    // Decide if this car will cause us to adjust speed
    td = computeTargetSpeed(last_s_processed, last_speed_processed, last_acc_processed, sd);
    trg_speed = td.final_speed;
    if (trg_speed - prev_trg_speed_ > 0 and trg_speed - prev_trg_speed_ < 3) {
      trg_speed = prev_trg_speed_;
    }
  } else {
    td = computeMatchTargetSpeed(last_s_processed, last_speed_processed, last_acc_processed, trg_speed);
  }

  double cur_speed, cur_acc, prev_ns;
  // Determine if the target speed has changed
  if (prev_trg_speed_ != trg_speed) {
    policy_      = MATCH_SPEED;
    cur_speed    = last_speed_processed;
    cur_acc      = last_acc_processed;
    prev_ns      = last_s_processed;
    //td = computeMatchTargetSpeed(last_s_processed, last_speed_processed, last_acc_processed, trg_speed);
    prev_trg_speed_ = trg_speed;
    computeTrajectory(td);
    if (DEBUG) {
      std::cout << "Computed new trajectory with " << forward_traj_.size() << " steps to target speed " << 
                   trg_speed << " in time " << td.time << " :\n"; 
      std::cout << " init_s: " << td.init_s << " init_speed: " << td.init_speed << 
                   " init_acc: " << td.init_acc << " final_s: " << td.final_s <<
                   " final_speed: " << td.final_speed << " final_acc: " << td.final_acc << std::endl;
    }
    // Clear stale trajectory data
    next_x_vals_.clear();
    next_y_vals_.clear();
    next_s_vals_.clear();
    prev_acc_.clear(); 
    prev_speed_.clear(); 
    sz = 0; // to make sure we flush the remaining unprocessed points
  // If target speed hasn't changed just continue with the previously computed trajectory
  } else {
    cur_acc   = prev_acc_.size()?    prev_acc_.back():    0;
    cur_speed = prev_speed_.size()?  prev_speed_.back():  0;
    prev_ns   = next_s_vals_.size()? next_s_vals_.back(): st.s;
    cur_traj_idx_ += processed;
  }

  double prev_speed = cur_speed;

  // Build the rest of the steps
  for (int i = sz; i < TRAJ_STEPS; ++i) {
    // If we completed our computed trajectory, just change policy to keep current speed
    if (policy_ != KEEP_SPEED and cur_traj_idx_ + i >= forward_traj_.size()) {
      policy_ = KEEP_SPEED; 
      if (DEBUG) std::cout << "POLICY: KEEP_SPEED  (TrajSize: " << forward_traj_.size() << ")\n";
    }

    double nx, ny, ns;
    if (policy_ == KEEP_SPEED) {
      ns = prev_ns + prev_speed * TIME_RES;  // TODO: perhaps use prev_speed?
      ns = fmod(ns, MAX_S);
      auto vecXY = getXYfromSpline(ns, 6);
      nx = vecXY[0];
      ny = vecXY[1];
    } else {
      nx = forward_traj_[cur_traj_idx_ + i].x; 
      ny = forward_traj_[cur_traj_idx_ + i].y; 
      ns = forward_traj_[cur_traj_idx_ + i].s;
    }
    next_x_vals_.push_back(nx);
    next_y_vals_.push_back(ny);
    next_s_vals_.push_back(ns);
    // Keep track of speed and acceleration
    if (ns < prev_ns) ns += MAX_S; // Wraparound case
    cur_speed  = (ns - prev_ns) / TIME_RES;
    cur_acc    = (cur_speed - prev_speed) / TIME_RES;
    prev_ns    = fmod(ns, MAX_S);
    prev_speed = cur_speed;
    prev_acc_.push_back(cur_acc);
    prev_speed_.push_back(cur_speed);
    if (DEBUG) {
      std::cout << "Added new trajectory step " << (i+1) << " with s: " << ns << " speed: " << cur_speed <<
                   " acc: " << cur_acc << std::endl;
    }   
  }

  if (false and DEBUG) {
    for (int i = 0; i < next_x_vals_.size(); ++i) {
      std::cout << "X: " << next_x_vals_[i] << " Y: " << next_y_vals_[i] << " S: " << next_s_vals_[i] << std::endl; 
    }
    std::cout << "=======================================================================\n";
    std::cout << "\n\n";
  }
  return;
}
