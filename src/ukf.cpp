#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //0.5 = 3.6km/h (sda_a * 2)
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // 0.39 = 2.8km/h
  std_yawdd_ = 0.39;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);  

  weights_ = VectorXd(2*n_aug_+1);

  P_ <<   0.1, 0.0, 0.0, 0.0, 0.0, 
          0.0, 0.1, 0.0, 0.0, 0.0, 
          0.0, 0.0, 0.5, 0.0, 0.0, 
          0.0, 0.0, 0.0, 1.5, 0.0, 
          0.0, 0.0, 0.0, 0.0, 0.6;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  //first measurement
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == meas_package.LASER) {
      x_ << meas_package.raw_measurements_[0], 
            meas_package.raw_measurements_[1], 
            4.0, 
            0.1, 
            0.01;
    } else if (meas_package.sensor_type_ == meas_package.RADAR) {
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      // Convert radar from polar to cartesian coordinates and initialize state.
      float px = rho * cos(phi);
      float py = rho * sin(phi);

      x_ << px, py, 4.0, 0.1, 0.01;
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }
  
  //convert from microsencod to second
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == meas_package.RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == meas_package.LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //the first column is reserved to mean
  Xsig_aug.col(0) = x_aug;

  //1) define sigma points

  //calculate square root of P_aug
  MatrixXd A = P_aug.llt().matrixL();

  for (int i = 0; i < n_aug_; i++) {
    MatrixXd sigmas = sqrt(lambda_ + n_aug_) * A.col(i);
    //the first half of the matrix is reserved to 'positives' directions    
    Xsig_aug.col(i+1) = x_aug + sigmas;
    //the second half of the matrix is reserved for directions opposite the first half
    Xsig_aug.col(i+1+n_aug_) = x_aug - sigmas;
  }

  //2) predict the sigma points

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //split the vectors in x and noise
    VectorXd x = Xsig_aug.col(i).head(5);
    VectorXd nu = Xsig_aug.col(i).tail(2);
    
    float v = x(2);
    float yaw = x(3);
    float yaw_rate = x(4);
    
    //noise vector
    
    float nu_a = nu(0);
    float nu_yaw_rate = nu(1);

    double px_p, py_p;
    //avoid division by zero    
    if (fabs(yaw_rate) > 0.0001) {
      px_p = v / yaw_rate * (sin(yaw + yaw_rate * delta_t) - sin(yaw));
      py_p = v / yaw_rate * (-cos(yaw + yaw_rate * delta_t) + cos(yaw)); 
    } else { //is zero    
      px_p = v * cos(yaw) * delta_t;
      py_p = v * sin(yaw) * delta_t; 
    }
    
    VectorXd pred = VectorXd(5);
    pred << px_p,
            py_p,
            0, 
            yaw_rate * delta_t, 
            0;
    
    VectorXd noise = VectorXd(5);
    noise << 0.5 * (delta_t * delta_t) * cos(yaw) * nu_a, 
             0.5 * (delta_t * delta_t) * sin(yaw) * nu_a, 
             delta_t * nu_a, 
             1.0/2 * (delta_t * delta_t) * nu_yaw_rate, 
             delta_t * nu_yaw_rate;
    
    Xsig_pred_.col(i) = x + pred + noise;    
  }

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //3) process model -  calculate the mean and covariance of the predicted state
  
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {  
  n_z_ = 2;

  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);

  //measurement model -  calculate the mean and covariance of the predicted state

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd state_vector = Xsig_pred_.col(i);
      float px = state_vector(0);
      float py = state_vector(1);
      
      Zsig.col(i) << px, py;
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //define a R measurement noise
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
                  
  //calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd diff = Zsig.col(i) - z_pred;
      S = S + weights_(i) * diff * diff.transpose();
  }
  S = S + R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred;
     
      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  VectorXd z = meas_package.raw_measurements_;

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  float nis_value = z_diff.transpose() * S.inverse() * z_diff;

  ofstream myfile;
  myfile.open("../data/nis_fusion.txt", fstream::out | fstream::app);

  myfile << nis_value << '\n';

  myfile.close();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  n_z_ = 3;
  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);

  //measurement model -  calculate the mean and covariance of the predicted state

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd state_vector = Xsig_pred_.col(i);
      float px = state_vector(0);
      float py = state_vector(1);
      float v = state_vector(2);
      float yaw = state_vector(3);
      
      float sqrt_pxpy = sqrt(px*px + py*py);
      
      Zsig.col(i) << sqrt_pxpy, 
                 atan2(py,px), 
                 (px*cos(yaw)*v+py*sin(yaw)*v) / sqrt_pxpy;
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //define a R measurement noise
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0, 
      0, 0, std_radrd_*std_radrd_;
                  
  //calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd diff = Zsig.col(i) - z_pred;

      //angle normalization
      while (diff(1)> M_PI) diff(1)-=2.*M_PI;
      while (diff(1)<-M_PI) diff(1)+=2.*M_PI;

      S = S + weights_(i) * diff * diff.transpose();
  }
  S = S + R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred;

      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  VectorXd z = meas_package.raw_measurements_;

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  float nis_value = z_diff.transpose() * S.inverse() * z_diff;

  ofstream myfile;
  myfile.open("../data/nis_fusion.txt", fstream::out | fstream::app);

  myfile << nis_value << '\n';

  myfile.close();
}
