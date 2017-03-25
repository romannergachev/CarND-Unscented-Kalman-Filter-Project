#include <iostream>
#include "ukf.h"

using namespace std;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

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

  is_initialized_ = false;

  n_aug_ = n_x_ + 2;

  //Xsig_pred_ = MatrixXd(n_aug_, n_aug_);

  lambda_ = 3 - n_aug_;

  x_aug_ = VectorXd(n_aug_);

  P_aug_ = MatrixXd(n_aug_, n_aug_);

  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);



  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);

  R_lidar_ = MatrixXd(n_z_lidar_, n_z_lidar_);

  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0,                       std_laspy_ * std_laspy_;

  P_ << 1, 0, 0,  0,  0,
        0, 1, 0,  0,  0,
        0, 0, 10, 0,  0,
        0, 0, 0,  10, 0,
        0, 0, 0,  0,  10;

  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /*****************************************************************************
 *  Initialization
 ****************************************************************************/
  if (!is_initialized_) {
    previous_timestamp_ = measurement_pack.timestamp_;
    /**
      * Initialize the state ekf_.x_ with the first measurement.
    */
    cout << "UKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double d_rho = measurement_pack.raw_measurements_[2];

      x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //extracting long intervals
  while (dt > 0.1) {
    Prediction(0.05);
    dt -= 0.05;
  }

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  //select correct n_z_ value
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    n_z_ = n_z_radar_;
  } else {
    n_z_ = n_z_lidar_;
  }

  z_pred_ = VectorXd(n_z_);

  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  S_ = MatrixXd(n_z_, n_z_);

  Tc_ = MatrixXd(n_x_, n_z_);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    try {
      PredictRadarMeasurement();
      UpdateState(measurement_pack);
      NIS_radar_ = ((measurement_pack.raw_measurements_-z_pred_).transpose())*
                    S_.inverse()*(measurement_pack.raw_measurements_-z_pred_);
    } catch (...) {

    }


  } else {
    // Laser updates
    PredictLidarMeasurement();
    UpdateState(measurement_pack);
    NIS_laser_ = ((measurement_pack.raw_measurements_-z_pred_).transpose())*
                  S_.inverse()*(measurement_pack.raw_measurements_-z_pred_);
  }


  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  SigmaPointPrediction(delta_t);
  PredictMeanAndCovariance();
}

void UKF::SigmaPointPrediction(double delta_t) {
  AugmentedSigmaPoints();
  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

void UKF::AugmentedSigmaPoints() {
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}


void UKF::PredictMeanAndCovariance() {
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff(3) = tools.normalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Predict measurement using a radar data
 */
void UKF::PredictRadarMeasurement() {
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    double sqRoot = sqrt(p_x * p_x + p_y * p_y);

    // measurement model
    Zsig_(0, i) = sqRoot;                        //r
    if (fabs(p_x) > 0.0001) {
      Zsig_(1, i) = atan(p_y/p_x);               //phi
    } else {
      throw 1;
    }


    if (fabs(sqRoot) > 0.0001) {
      Zsig_(2, i) = (p_x * v1 + p_y * v2) / sqRoot;
    } else {
      throw 1;
    }
  }

  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    z_diff(1) = tools.normalizeAngle(z_diff(1));

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  S_ += R_radar_;
}

/**
 * Predict measurement using a laser data
 */
void UKF::PredictLidarMeasurement() {
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // measurement model
    Zsig_(0, i) = Xsig_pred_(0, i);
    Zsig_(1, i) = Xsig_pred_(1, i);
  }

  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  S_ += R_lidar_;
}

/**
  * Updates the state and the state covariance matrix
  */
void UKF::UpdateState(MeasurementPackage meas_package) {
  //calculate cross correlation matrix
  Tc_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    z_diff(1) = tools.normalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff(3) = tools.normalizeAngle(x_diff(3));

    Tc_ += weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc_ * S_.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred_;

  //angle normalization
  z_diff(1) = tools.normalizeAngle(z_diff(1));



  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_ * K.transpose();
}
