#include <cmath>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::Map;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = Matrix<double, 5, 1>();

  // initial covariance matrix
  P_ = Matrix<double, 5, 5>::Identity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  lambda_ = -2.0;
  weights_ = Matrix<double, n_pts_, 1>::Constant(1.0 / 2.0 / (lambda_ + n_aug_));
  weights_(0) *= 2*lambda_;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {

  const VectorXd& measurements = meas_package.raw_measurements_;
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << measurements[0],
           measurements[1],
           0.0,
           0.0,
           0.0;
    } else { // Initialize with Radar measurement.
      float x_comp = std::cos(measurements[1]);
      float y_comp = std::sin(measurements[1]);
      x_ << measurements[0] * x_comp,
            measurements[0] * y_comp,
            0.0,
            0.0,
            0.0;
    }
    is_initialized_ = true;
  }
  else {
    UKF::Prediction(meas_package.timestamp_ - previous_timestamp_);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UKF::UpdateLidar(measurements);
    } else {
        UKF::UpdateRadar(meas_package);
    }
  }
  previous_timestamp_ = meas_package.timestamp_;
}

void UKF::Prediction(double delta_t) {
  // Create augmented covariance matrix and x.
  Matrix<double, n_aug_, n_aug_> P_aug = Matrix<double, n_aug_, n_aug_>::Zero();
  P_aug.topLeftCorner<5, 5>() = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_; 

  Matrix<double, n_aug_, 1> x_aug;
  x_aug << x_, 0.0, 0,0;

  // Calculate sigma points.
  double coef = std::sqrt(lambda_ + n_aug_);
  Matrix<double, n_aug_, n_aug_> sigmas = 
    coef * Matrix<double, n_aug_, n_aug_>(P_aug.llt().matrixL());
  Matrix<double, n_aug_, n_pts_> Xsig_aug;
  Xsig_aug << x_aug, sigmas.colwise() + x_aug, (-1 * sigmas).colwise() + x_aug;

  // Calculate predictions for the sigma points.
  double delta_t2 = delta_t * delta_t;

  for (int i = 0; i < n_pts_; ++i) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yaw_rate = Xsig_aug(4,i);
    double noise_a = Xsig_aug(5,i);
    double noise_ay = Xsig_aug(6,i);
        
    // Precompute some reused terms
    double cos_yaw = cos(yaw);
    double sin_yaw = sin(yaw);
    
    // Initialize with previous value and noise.
    double new_px = px + 0.5*delta_t2*cos_yaw*noise_a;
    double new_py = py + 0.5*delta_t2*sin_yaw*noise_a;
    double new_v = v + delta_t * noise_a;
    double new_yaw = yaw + 0.5*delta_t2*noise_ay;
    double new_yaw_rate = yaw_rate + delta_t * noise_ay;
    // Add linear-model or curved model terms
    if (yaw_rate == 0.0) {
      new_px += v*cos_yaw*delta_t;
      new_py += v*sin_yaw*delta_t;
            
    } else {
      new_px += v/yaw_rate*(sin(yaw + yaw_rate*delta_t) - sin_yaw);
      new_py += v/yaw_rate*(-1.0*cos(yaw + yaw_rate*delta_t) + cos_yaw);
      new_yaw += yaw_rate * delta_t;
    }
    Xsig_pred_(0,i) = new_px;
    Xsig_pred_(1,i) = new_py;
    Xsig_pred_(2,i) = new_v;
    Xsig_pred_(3,i) = new_yaw;
    Xsig_pred_(4,i) = new_yaw_rate;

    // TODO Assertion failing here!
    x_ = (Xsig_pred_.array().rowwise() *
          weights_.transpose().array()).rowwise().sum();
    Matrix<double, n_x_, n_x_> P = Matrix<double, n_x_, n_x_>::Zero();
    for (int i = 0; i < n_pts_; ++i) {
      Matrix<double, n_x_, 1> residual = Xsig_pred_.col(i) - x_;
      P += weights_(i) * residual * residual.transpose();
    }
    P_ = P;
  }
}

void UKF::UpdateLidar(const VectorXd&  measurements) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  const Matrix<double, 2, n_pts_> Zsig = 
       Xsig_pred_.block<2, n_pts_>(0,0);
   Vector2d z_pred = Zsig * weights_;
  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}
