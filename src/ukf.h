#ifndef UKF_H
#define UKF_H

#include <fstream>

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage& meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const Eigen::VectorXd& measurements);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const Eigen::VectorXd& measurements);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::Matrix<double, 5, 1> x_;

  // state covariance matrix
  Eigen::Matrix<double, 5, 5> P_;

  // time when the state is true, in us
  long long time_us_;

  long long previous_timestamp_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Directory location to write NIS data.
  std::ofstream nis_file_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;
  double var_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;
  double var_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;
  double var_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;
  double var_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;
  double var_radrd_ ;

  // State dimension
  static constexpr int n_x_ = 5;

  // Augmented state dimension
  static constexpr int n_aug_ = n_x_ + 2;

  // Number of sigma points.
  static constexpr int n_pts_ = 2*n_aug_ + 1;

  // predicted sigma points matrix
  Eigen::Matrix<double, 5, n_pts_> Xsig_pred_;

  // Weights of sigma points
  Eigen::Matrix<double, n_pts_, 1>  weights_;

  // Sigma point spreading parameter
  double lambda_;
};

#endif  // UKF_H
