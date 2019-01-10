#include <cmath>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::Matrix;
using Eigen::MatrixXd;
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
}

void UKF::Prediction(double delta_t) {
  // Create augmented covariance matrix and x.
  constexpr int n = 5;
  constexpr int n_aug = n+2;
  Matrix<double, n_aug, n_aug> P_aug = Matrix<double, n_aug, n_aug>::Zero();
  P_aug.topLeftCorner<n, n>() = P_;
  P_aug(n,n) = std_a_ * std_a_;
  P_aug(n+1,n+1) = std_yawdd_ * std_yawdd_; 

  Matrix<double, n_aug, 1> x_aug;
  x_aug << x_, 0.0, 0,0;

  // Calculate sigma points.
  double coef = std::sqrt(lambda_ + n_aug);
  Matrix<double, n_aug, n_aug> sigmas = 
    coef * Matrix<double, n_aug, n_aug>(P_aug.llt().matrixL());
  constexpr int n_pts = 2*n_aug + 1;
  Matrix<double, n_aug, n_pts> Xsig_aug;
  Xsig_aug << x_aug, sigmas.colwise() + x_aug, (-1 * sigmas).colwise() + x_aug;

  // Calculate predictions for the sigma points.
  Matrix<double, n, n_pts> Xsig_pred;
  double delta_t2 = delta_t * delta_t;

  for (int i = 0; i < n_pts; ++i) {
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
    Xsig_pred(0,i) = new_px;
    Xsig_pred(1,i) = new_py;
    Xsig_pred(2,i) = new_v;
    Xsig_pred(3,i) = new_yaw;
    Xsig_pred(4,i) = new_yaw_rate;

    Matrix<double, n_pts, 1>  weights = 
       Matrix<double, n_pts, 1>::Constant(1.0 / 2.0 / (lambda_ + n_aug));
    weights(0) *= 2*lambda_;
    x_ = (Xsig_pred.array().rowwise() *
          weights.transpose().array()).rowwise().sum();
    Matrix<double, n, n> P = Matrix<double, n, n>::Zero();
    for (int i = 0; i < n_pts; ++i) {
      Matrix<double, n, 1> residual = Xsig_pred.col(i) - x_;
      P += weights(i) * residual * residual.transpose();
    }
    P_ = P;
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}
