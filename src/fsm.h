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

static constexpr bool   DEBUG      = true;
static constexpr double MPH2MPS    = 0.44704;  // MPH to m/s
static constexpr double TIME_RES   = 0.02;     // Time resolution in sec
static constexpr int    TRAJ_STEPS = 50;       // trajectory steps
static constexpr double MAX_SPEED  = 50 * MPH2MPS; // in m/sec
static constexpr double MAX_ACC    = 10;       // in m/sec^2 
static constexpr double MAX_JERK   = 10;       // Max jerk in m/sec^3  
static constexpr double MAX_S      = 6945.554; // Max S (length) of our track

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
  double acc;
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
  double init_speed;
  double init_acc;
  double final_s;
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

class FSM {
private:
  typedef enum { KEEP_LANE, MOVE_LEFT, MOVE_RIGHT } Policy;
  Policy                policy_;
  std::vector<double>   prev_acc_;
  std::vector<double>   prev_speed_;
  std::vector<double>   next_x_vals_;
  std::vector<double>   next_y_vals_;
  std::vector<double>   next_s_vals_;
  const Map &           map_;
  SplineData            spd_;
  double                prev_trg_speed_;
  size_t                cur_traj_idx_;
  std::vector<TrajNode> forward_traj_;

  // How much delta_speed will we shed when decelerating from MAX_ACC to 0?
  double                delta_speed_to_zero_acc_;

  void generateFullSplines();
  std::vector<double> getXYfromSpline(double s, double d_offset);
  TrajData computeAccelerateToTrgSpeed(double init_s, double init_speed, double init_acc, double trg_speed);
  void computeTrajectory(TrajData & traj_data);

public:
  FSM(const Map & map);
  void update(State & st);
  std::pair<std::vector<double>, std::vector<double>> getNextXYVals();
};

#endif  // FSM_H
