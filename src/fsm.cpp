#include "fsm.h"
#include <iostream>
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
  policy_(KEEP_LANE), prev_acc_(TRAJ_STEPS, 0), prev_speed_(TRAJ_STEPS, 0), map_(map), prev_trg_speed_(0),
  cur_traj_idx_(0), delta_speed_to_zero_acc_(0)
{
  generateFullSplines();

  // Compute some const information
  int shed_acc_steps = ceil(MAX_ACC / MAX_JERK / TIME_RES);
  double a_shed = 0;
  for (int i = 0; i < shed_acc_steps; ++i) {
    a_shed += MAX_JERK * TIME_RES;
    delta_speed_to_zero_acc_ += a_shed * TIME_RES;
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

// Return copies of internal next XY values
std::pair<vector<double>, vector<double>> FSM::getNextXYVals()
{
  return {next_x_vals_, next_y_vals_};
} 


// Method to compute and return the minimum amount of road length and time
// it would take to acceleratate to target speed with 0 final acceleration
TrajData FSM::computeAccelerateToTrgSpeed(double init_s, double init_speed, double init_acc, double trg_speed) 
{
  // How many timesteps do we need at least to shed our initial acceleration?
  int shed_acc_steps      = ceil(init_acc / MAX_JERK / TIME_RES);
  double delta_speed      = trg_speed - init_speed;
  double delta_speed_shed = 0;
  double a_shed           = 0;
  // How much speed are we going to add while shedding initial acceleration?
  for (int i = 0; i < shed_acc_steps; ++i) {
    a_shed += MAX_JERK * TIME_RES;
    delta_speed_shed += a_shed * TIME_RES;
  }
  // Update the desired delta_speed
  delta_speed -= delta_speed_shed;
  // The flip point, to go from accelerating to decelerating
  // is the new delta_speed / 2
  double flip_point = (delta_speed > 0)? delta_speed / 2: 0;
  double T = 0;
  double S = init_s;
  double cur_speed = init_speed;
  double cur_acc = init_acc;
  while (cur_speed < trg_speed) {
    std::cout << "Cur speed: " << cur_speed << " Cur acc: " << cur_acc << " Cur S: " << S <<
                 " Flip point: " << flip_point << std::endl;
    T          += TIME_RES;
    S          += cur_speed * TIME_RES;
    cur_speed  += cur_acc * TIME_RES;
    if (cur_speed > MAX_SPEED) cur_speed = MAX_SPEED;
    double jerk = ((trg_speed - cur_speed) > flip_point)? MAX_JERK: -MAX_JERK;
    cur_acc    += jerk * TIME_RES;
    if (cur_acc > MAX_ACC) {
      cur_acc = MAX_ACC;
      // Update flip point as well
      flip_point = delta_speed_to_zero_acc_;
    }
    if (cur_acc < 0) {
      cur_acc = 0;
      break;
    }
  }

  TrajData ret;
  ret.init_s      = init_s;
  ret.init_speed  = init_speed;
  ret.init_acc    = init_acc;
  ret.final_s     = S;
  ret.final_speed = trg_speed;
  ret.final_acc   = cur_acc;
  ret.time        = T;

  return ret;
}

void FSM::computeTrajectory(TrajData & traj_data)
{
  forward_traj_.clear();
  cur_traj_idx_ = 0;
  auto coeffs = computePolyCoefficients(traj_data);
  // Generate trajector until we reach the required time
  double time = 0;
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

void FSM::update(State & st) 
{
  int sz = st.path_x.size();

  // How many points left in the previous path?
  if (sz == TRAJ_STEPS) {
    if (DEBUG) std::cout << "Returning same lists since nothing was processed\n";
    return; // nothing processed, the lists remain the same
  }

  next_x_vals_.clear();
  next_y_vals_.clear();

  // First copy over the old trajectory
  if (sz) {
    next_x_vals_.insert(next_x_vals_.begin(), st.path_x.begin(), st.path_x.end());
    next_y_vals_.insert(next_y_vals_.begin(), st.path_y.begin(), st.path_y.end());
  }

  // TODO: determine target speed based on sensor fusion
  double trg_speed = MAX_SPEED;

  int processed = TRAJ_STEPS - sz;
  if (DEBUG) {
    std::cout << "=======================================================================\n";
    std::cout << "Processed: " << (sz? std::to_string(processed): "N/A") << std::endl;
  }
  if (sz) next_s_vals_.erase(next_s_vals_.begin(), next_s_vals_.begin() + processed);

  // Erase processed first elements 
  prev_acc_.erase(prev_acc_.begin(), prev_acc_.begin() + processed);
  prev_speed_.erase(prev_speed_.begin(), prev_speed_.begin() + processed);
  double latest_acc = prev_acc_.size()? prev_acc_.back(): 0;
  double latest_speed = prev_speed_.size()? prev_speed_.back(): 0;
  double first_acc = prev_acc_.size()? prev_acc_.front(): 0;
  double first_speed = prev_speed_.size()? prev_speed_.front(): 0;

  if (DEBUG) {
    std::cout << "Current state (s,x,y,speed,acc,path sz): (" << st.s << ", " << st.x << ", " << st.y <<
                 ", " << first_speed << ", " << first_acc << ", " << sz << ")\n";
  }

  double cur_s      = (sz)? st.end_s: st.s;

  // Determine if the target speed has changed
  if (prev_trg_speed_ != trg_speed) {
    auto td = computeAccelerateToTrgSpeed(cur_s, first_acc, first_speed, trg_speed);
    if (DEBUG) {
      std::cout << "Computed new trajectory to target speed " << trg_speed << " in time " << td.time << " :\n"; 
      std::cout << " init_s: " << td.init_s << " init_speed: " << td.init_speed << 
                   " init_acc: " << td.init_acc << " final_s: " << td.final_s <<
                   " final_speed: " << td.final_speed << " final_acc: " << td.final_acc << std::endl;
    }
    prev_trg_speed_ = trg_speed;
    computeTrajectory(td);
    // Clear stale trajectory data
    next_x_vals_.clear();
    next_y_vals_.clear();
    next_s_vals_.clear();
    prev_acc_.clear(); 
    prev_speed_.clear(); 
    sz = 0; // to make sure we flush the remaining unprocessed points
  } else {
    cur_traj_idx_ += processed;
  }

  double cur_speed  = latest_speed;
  double cur_acc    = latest_acc;
  double prev_ns    = cur_s;
  double prev_speed = cur_speed;

  // Build the rest of the steps
  for (int i = sz; i < TRAJ_STEPS; ++i) {
    // If we completed our computed trajectory, just keep current speed
    const bool keep_speed = (cur_traj_idx_ + i) >= forward_traj_.size(); 
    double nx, ny, ns;
    if (keep_speed) {
      ns = prev_ns + prev_trg_speed_ * TIME_RES;  // TODO: perhaps use prev_speed?
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
  }

  if (DEBUG) {
    for (int i = 0; i < next_x_vals_.size(); ++i) {
      std::cout << "X: " << next_x_vals_[i] << " Y: " << next_y_vals_[i] << " S: " << next_s_vals_[i] << std::endl; 
    }
    std::cout << "=======================================================================\n";
    std::cout << "\n\n";
  }
  return;
}
