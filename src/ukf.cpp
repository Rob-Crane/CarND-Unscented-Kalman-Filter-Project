#include <cmath>
#include <cstdlib>
#include <string>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::Map;
using Eigen::Matrix;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Vector2d;
using Eigen::Vector3d;
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
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // If NIS_DATA_DIR environment variable set, write NIS to data file.
  if (const char* nis_dir = std::getenv("NIS_DATA_DIR")) {
    char buffer[30];
    std::sprintf(buffer, "a%.2f_%.2f.csv", std_a_, std_yawdd_);
    std::string nis_filename = std::string(nis_dir) + "/" + std::string(buffer);
    nis_file_.open(nis_filename);
  }

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
  var_laspx_ = std_laspx_ * std_laspx_;
  var_laspy_ = std_laspy_ * std_laspy_;

  var_radr_ = std_radr_ * std_radr_;
  var_radphi_ = std_radphi_ * std_radphi_;
  var_radrd_ = std_radrd_ * std_radrd_;
}


UKF::~UKF() {
  nis_file_.close();
}


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
    double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1.0E6;
    UKF::Prediction(delta_t);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UKF::UpdateLidar(measurements);
    } else {
        UKF::UpdateRadar(measurements);
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
  x_aug << x_, 0.0, 0.0;

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
  }
  x_ = (Xsig_pred_.array().rowwise() *
        weights_.transpose().array()).rowwise().sum();
  Matrix<double, n_x_, n_x_> P = Matrix<double, n_x_, n_x_>::Zero();
  for (int i = 0; i < n_pts_; ++i) {
    Matrix<double, n_x_, 1> residual = Xsig_pred_.col(i) - x_;
    P += weights_(i) * residual * residual.transpose();
  }
  P_ = P;
}

void UKF::UpdateLidar(const VectorXd&  measurements) {
  // TODO: calculate the lidar NIS
  const Matrix<double, 2, n_pts_> Zsig = 
       Xsig_pred_.block<2, n_pts_>(0,0);
   Vector2d z_pred = Zsig * weights_;
   Matrix<double, 2, n_pts_> Zres = Zsig.colwise() - z_pred;
   Matrix2d S = Matrix2d::Zero();
   for (int i = 0; i < n_pts_; ++i) {
     Vector2d resid = Zres.col(i);
     S+= weights_(i) * resid * resid.transpose();
   }
   S(0,0) += var_laspx_;
   S(1,1) += var_laspy_;

   Matrix<double, n_x_, n_pts_> Xres = Xsig_pred_.colwise() - x_;
   Matrix<double, n_x_, 2> T = Xres * weights_.asDiagonal() * Zres.transpose();
   Matrix<double, n_x_, 2> K = T * S.inverse();

   Vector2d innovation = measurements - z_pred;
   // Calculate normalized innovation squared (NIS)
   double nis = innovation.transpose() * S.inverse() * innovation;
   if (nis_file_.is_open()) {
     nis_file_ << "LASER," << std::to_string(nis) << "," << std::endl;
   }

   x_ += K*(innovation);
   P_ -= K*S*K.transpose();
}

void UKF::UpdateRadar(const VectorXd&  measurements) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  Matrix<double, 3, n_pts_> Zsig;
  for (int i = 0; i < n_pts_; ++i) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    Zsig(0, i) = sqrt(px*px + py*py);
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = v*(px*cos(yaw) + py*sin(yaw)) / Zsig(0, i);
  }
  Vector3d z_pred = Zsig*weights_;
  Matrix<double, 3, n_pts_> Zres = Zsig.colwise() - z_pred;
  Matrix3d S = Matrix3d::Zero();
  for (int i = 0; i < n_pts_; ++i) {
    Vector3d resid = Zres.col(i);
    S+= weights_(i) * resid * resid.transpose();
  }
  S(0,0) += var_radr_; 
  S(1,1) += var_radphi_; 
  S(2,2) += var_radrd_; 
  Matrix<double, n_x_, n_pts_> Xres = Xsig_pred_.colwise() - x_;
  Matrix<double, n_x_, 3> T = Xres * weights_.asDiagonal() * Zres.transpose();
  Matrix<double, n_x_, 3> K = T * S.inverse();
  Vector3d innovation = measurements - z_pred;
  // Calculate normalized innovation squared (NIS)
  double nis = innovation.transpose() * S.inverse() * innovation;
   if (nis_file_.is_open()) {
     nis_file_ << "RADAR," << std::to_string(nis) << "," << std::endl;
   }
  x_ += K*innovation;
  P_ -= K*S*K.transpose();
}
