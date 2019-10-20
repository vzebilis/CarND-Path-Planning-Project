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
static constexpr double LANE_CHANGE_MAR = 7;        // Absolute margin for lane change
static constexpr double LANE_CHANGE_FCT = LANE_WIDTH / 0.1; // Factor for computing the time to change lanes
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
} Context;

typedef struct {
  double x;
  double y;
  double s;
  double d;
  double v;
  double a;
} PositionData;

typedef struct {
  tk::spline x;
  tk::spline y;
  tk::spline dx;
  tk::spline dy;
  tk::spline s;
  tk::spline d;
} SplineData;

typedef struct {
  double init_s;
  double init_d;
  double init_speed;
  double init_acc;
  double mid_s;
  double mid_d;
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

typedef enum { KEEP_LANE, MATCH_SPEED, CHANGE_LANE } Policy;
typedef enum { CENTER, RIGHT, LEFT } Lane;

typedef struct {
  Policy pol;
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
    PositionData          p_;
    double                trg_speed_;
    size_t                cur_traj_idx_;
    std::vector<TrajNode> forward_traj_;
public:
    State(FSM & fsm, double trg_speed) : fsm_(fsm), trg_speed_(trg_speed), cur_traj_idx_(0) {
      trg_speed_ = fmin(trg_speed_, MAX_SPEED * MAX_SPEED_MAR);
    }
    virtual ~State() {}
    virtual void update() = 0;
    virtual Policy nextPolicy() = 0;
    virtual TrajData computeTargetPos(PositionData & p, SensorData * sd) { return {}; }
    virtual void computeTrajectory(TrajData & td);
    // NON-VIRTUAL
    double getTargetSpeed() const { return trg_speed_; }
};

class KeepLane : public State {
public:
    using State::State; // using base ctor
    void update() override;
    Policy nextPolicy() override { return KEEP_LANE; }
};

class ChangeLane : public State {
public:
    ChangeLane(FSM & fsm, double trg_speed, SensorData * sd);
    void update() override;
    Policy nextPolicy() override;
    TrajData computeTargetPos(PositionData & p, SensorData * sd) override;
    void computeTrajectory(TrajData & td) override;
};

class MatchSpeed : public State {
public:
    MatchSpeed(FSM & fsm, double trg_speed, SensorData * sd);
    void update() override;
    Policy nextPolicy() override;
    TrajData computeTargetPos(PositionData & p, SensorData * sd) override;
};

class FSM {
public:
  State *               state_;
  bool                  first_time_;
  Context *             cxt_;
  Policy                prev_pol_;
  Lane                  lane_;

  const Map &           map_;
  SplineData            spd_;

  std::vector<double>   prev_acc_;
  std::vector<double>   prev_speed_;
  std::vector<double>   next_x_vals_;
  std::vector<double>   next_y_vals_;
  std::vector<double>   next_s_vals_;
  std::vector<double>   next_d_vals_;
  double                prev_trg_speed_;
  // How much delta_speed will we shed when decelerating from MAX_ACC to 0?
  double                delta_speed_to_zero_acc_;

  void generateFullSplines();
  SplineData * getShortSplines(TrajData & td);
  std::vector<double> getXYfromSpline(double s, double d_offset);
  TrajData computeMatchTargetSpeed(double init_s, double init_d, double init_speed, double init_acc,
                                   double final_d, double trg_speed);
  bool canChangeLane(const PositionData & p, SensorData & ahead, SensorData & behind);
  SensorData checkLanes(const PositionData & p);
  void updateLane(double d);
  void handleFirstTime();

public:
  FSM(const Map & map);
  void update(Context & cxt);
  void clearProcessed();
  void clearAll() {
      prev_acc_.clear(); prev_speed_.clear();
      next_x_vals_.clear(); next_y_vals_.clear();
      next_s_vals_.clear(); next_d_vals_.clear();
  }
  PositionData getLastProcessed() const;
  PositionData getMostRecentNotProcessed() const;
  void incrementByProcessed(size_t & val) const;
  void applyTrajectory(PositionData & p, std::vector<TrajNode> & forward_traj, size_t cur_traj_idx);
  std::pair<std::vector<double>, std::vector<double>> getNextXYVals() { return {next_x_vals_, next_y_vals_}; }
};

#endif  // FSM_H
