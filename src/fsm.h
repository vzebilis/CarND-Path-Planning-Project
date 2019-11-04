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
static constexpr int    TRAJ_STEPS      = 50;       // trajectory steps to generate each time
static constexpr double MAX_SPEED       = 50 * MPH2MPS; // specified in MPH, converted to m/sec
static constexpr double MAX_ACC         = 0.95 * 10;// in m/sec^2
static constexpr double MAX_JERK        = 10;       // Max jerk in m/sec^3
static constexpr double MAX_S           = 6945.554; // Max S (length) of our track
static constexpr double MAX_SPEED_MAR   = 0.98;     // MAX_SPEED margin factor
static constexpr double FRONT_CAR_MAR   = 10;       // Safety margin in meters from the front car
static constexpr double LANE_CHANGE_MAR = 7;        // Absolute margin in meters for lane change
static constexpr double SENSOR_MAX_DIST = 30;       // Maximum distance to look ahead for sensor data

static constexpr double LANE_WIDTH      = 4;        // Lane width in meters
static constexpr double LANE_CHANGE_FCT = 40;       // Factor for computing the time to do lane changing
static constexpr double HALF_LANE_WIDTH = LANE_WIDTH / 2;
static constexpr double LEFT_LANE_D     = HALF_LANE_WIDTH;
static constexpr double CENTER_LANE_D   = HALF_LANE_WIDTH + LANE_WIDTH;
static constexpr double RIGHT_LANE_D    = HALF_LANE_WIDTH + 2 * LANE_WIDTH;
static constexpr double STARTING_LANE   = CENTER_LANE_D; // Starting lane

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
} Context;

typedef struct {
  double x;
  double y;
  double s;
  double d;
} PositionData;

typedef struct {
  double yaw;
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} Movement;

typedef struct {
  tk::spline x;
  tk::spline y;
  tk::spline dx;
  tk::spline dy;
  tk::spline s;
  tk::spline d;
} SplineData;

typedef struct {
  double x;
  double y;
} TrajNode;

typedef enum { KEEP_LANE, MATCH_SPEED, CHANGE_LANE } Policy;
typedef enum { CENTER, RIGHT, LEFT } Lane;

typedef struct {
  Policy policy;
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

class FSM;

class State {
protected:
    FSM &                 fsm_;
    double                trg_speed_;
    size_t                cur_traj_idx_;
    std::vector<TrajNode> forward_traj_;
public:
    State(FSM & fsm, double trg_speed) : fsm_(fsm), trg_speed_(trg_speed), cur_traj_idx_(0) {
      trg_speed_ = fmin(trg_speed_, MAX_SPEED * MAX_SPEED_MAR);
    }
    virtual ~State() {}
    virtual Policy nextPolicy() = 0;
    // NON-VIRTUAL
    void update();
    double getTargetSpeed() const { return trg_speed_; }
};

class KeepLane : public State {
public:
    using State::State; // using base ctor
    Policy nextPolicy() override { return KEEP_LANE; }
};

class ChangeLane : public State {
public:
    ChangeLane(FSM & fsm, double trg_speed, SensorData * sd);
    Policy nextPolicy() override;
};

class MatchSpeed : public State {
public:
    MatchSpeed(FSM & fsm, double trg_speed, SensorData * sd);
    Policy nextPolicy() override;
    void matchTargetSpeed(double trg_v);
};

class FSM {
public:
  State *               state_;
  Context *             cxt_;
  Policy                cur_policy_;
  Lane                  lane_;

  const Map &           map_;

  std::vector<double>   next_x_vals_;
  std::vector<double>   next_y_vals_;
  std::vector<double>   next_s_vals_;
  std::vector<double>   next_d_vals_;
  double                cur_trg_speed_;
  bool                  first_time_;

  tk::spline * createLocalSpline(double d, const Movement & prev_m, double & ref_x, double & ref_y, double & ref_yaw);
  void createXYFromSpline(double d, std::vector<TrajNode> & ft, size_t & cti);
  bool canChangeLane(SensorData & ahead, SensorData & behind);
  SensorData checkLanes(const PositionData & p);
  void updateLane(double d);
  void handleFirstTime();

public:
  FSM(const Map & map);
  void update(Context & cxt);
  void clearProcessed();
  void clearAll() {
      next_x_vals_.clear(); next_y_vals_.clear();
      next_s_vals_.clear(); next_d_vals_.clear();
  }
  Movement getLastPlannedMovement();
  PositionData getLastProcessed() const;
  PositionData getMostRecentNotProcessed() const;
  void incrementByProcessed(size_t & val) const;
  void applyTrajectory(PositionData & p, std::vector<TrajNode> & forward_traj, size_t & cur_traj_idx);
  std::pair<std::vector<double>, std::vector<double>> getNextXYVals() { return {next_x_vals_, next_y_vals_}; }
};

#endif  // FSM_H
