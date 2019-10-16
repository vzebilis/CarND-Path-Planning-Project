#ifndef FSM_H
#define FSM_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
// Using recommended Cubic Spline implementation from:
// https://kluge.in-chemnitz.de/opensource/spline/
#include "spline.h"

static constexpr bool   DEBUG           = true;
static constexpr double MPH2MPS         = 0.44704;  // MPH to m/s
static constexpr double TIME_RES        = 0.02;     // Time resolution in sec
static constexpr int    TRAJ_STEPS      = 50;       // trajectory steps
static constexpr double MAX_SPEED       = 50 * MPH2MPS; // specified in MPH, converted to m/sec
static constexpr double MAX_ACC         = 10;       // in m/sec^2 
static constexpr double MAX_JERK        = 10;       // Max jerk in m/sec^3  
static constexpr double MAX_S           = 6945.554; // Max S (length) of our track
static constexpr double MAX_SPEED_MAR   = 0.99;     // MAX_SPEED margin factor 
static constexpr double LANE_WIDTH      = 4;        // Lane width in meters
static constexpr double FRONT_CAR_MAR   = 10;       // Safety margin from the front car (m)
static constexpr double SENSOR_MAX_DIST = 30;       // Maxmum distance to look ahead for sensor data
static constexpr double STARTING_LANE   = 6;        // Starting lane. 2: left, 6: middle, 10: right

typedef struct {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> s;
  std::vector<double> dx;
  std::vector<double> dy;
} Map;

typedef struct {
  double x;
  double y;
  double s;
  double d;
  double yaw;
  double speed;
  nlohmann::basic_json<>::value_type path_x;
  nlohmann::basic_json<>::value_type path_y;
  double end_s;
  double end_d;
  nlohmann::basic_json<>::value_type sensor_fusion;
} State;

typedef struct {
  tk::spline sx;
  tk::spline sy;
  tk::spline sdx;
  tk::spline sdy;
} SplineData;

typedef struct {
  double init_s;
  double init_d;
  double init_speed;
  double init_acc;
  double final_s;
  double final_d;
  double final_speed;
  double final_acc;
  double time;
} TrajData;

typedef struct {
  double x;
  double y;
  double s;
  double d;
} TrajNode;

typedef struct {
  int    id;
  double x;
  double y;
  double vx;
  double vy;
  double v;
  double s;
  double d;
  double diff_s;
} SensorData;

class FSM {
private:
  typedef enum { KEEP_LANE, MOVE_LEFT, MOVE_RIGHT } Policy;
  typedef enum { CENTER, RIGHT, LEFT } Lane;
  Policy                policy_;
  Policy                prev_policy_;
  Lane                  lane_;
  std::vector<double>   prev_acc_;
  std::vector<double>   prev_speed_;
  std::vector<double>   next_x_vals_;
  std::vector<double>   next_y_vals_;
  std::vector<double>   next_s_vals_;
  std::vector<double>   next_d_vals_;
  const Map &           map_;
  SplineData            spd_;
  double                prev_trg_speed_;
  size_t                cur_traj_idx_;
  std::vector<TrajNode> forward_traj_;

  // How much delta_speed will we shed when decelerating from MAX_ACC to 0?
  double                delta_speed_to_zero_acc_;

  void generateFullSplines();
  std::vector<double> getXYfromSpline(double s, double d_offset);
  TrajData computeMatchTargetSpeed(double init_s, double init_d, double init_speed, double init_acc, 
                                   double final_d, double trg_speed);
  //TrajData computeAccelerateToTrgSpeed(double init_s, double init_speed, double init_acc, double trg_speed);
  void computeTrajectory(TrajData & traj_data);
  TrajData computeTargetSpeed(double init_s, double init_d, double init_speed, double init_acc, SensorData & sd);
  SensorData checkLanes(State & st);
  void updateLane(double d);

public:
  FSM(const Map & map);
  void update(State & st);
  std::pair<std::vector<double>, std::vector<double>> getNextXYVals() { return {next_x_vals_, next_y_vals_}; }
};

#endif  // FSM_H
