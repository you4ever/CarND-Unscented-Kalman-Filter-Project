#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#define EPS 0.00001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // cout << "initialize UKF" << endl;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  is_initialized_ = false;

  // size of state vector
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.12;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.12;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.05;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties. 

  Hint: one or more values initialized above might be wildly off...
  */

  // Size of augmented state vector
  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_aug_;

  // cout << "   Done initialize UKF" << endl;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if (!(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) &&
		!(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {
		return;
	}

	// cout << "ProcessMeasurement" << endl;
	if (!is_initialized_) {

		// Initialize covairance matrix of state vector
		P_ << 1, 0, 0, 0, 0,
			  0, 1, 0, 0, 0,
			  0, 0, 2, 0, 0,
			  0, 0, 0, M_PI/2, 0,
			  0, 0, 0, 0, M_PI/2;

		// set weights for average sigma points
		// create vector for weights
		weights_ = VectorXd(2 * n_aug_ + 1); 
		weights_.fill(0.5 / (lambda_ + n_aug_));
		weights_(0) = lambda_ / (lambda_ + n_aug_);

		// Initialize measurements
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = meas_package.raw_measurements_[0]; // range
			float phi = meas_package.raw_measurements_[1]; // bearing
			float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho
																   // Coordinates convertion from polar to cartesian
			float px = rho * cos(phi);
			float py = rho * sin(phi);
			float vx = rho_dot * cos(phi);
			float vy = rho_dot * sin(phi);
			//float v = sqrt(vx * vx + vy * vy);

			x_ << px, py, rho_dot, 0., 0. ;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			float px = meas_package.raw_measurements_[0];
			float py = meas_package.raw_measurements_[1];
			x_ << px, py, 0., 0., 0.;
		}
		// Deal with the special case initialisation problems

		if (fabs(x_(0)) < EPS && fabs(x_(1)) < EPS){
			x_(0) = EPS;
			x_(1) = EPS;
		}
		// Save the initiall timestamp for dt calculation
		time_us_ = meas_package.timestamp_;

		Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
		Xsig_pred_.fill(0.0);

		// Done initializing, no need to predict or update
		is_initialized_ = true;

		return;
	}

	//double timestamp_us = meas_package.timestamp_ / 1e6;
	// Calculate time elapsed from last measurement in seconds
	double dt_s = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;

	// Reset filter if too long has elpsed since last measurement
	if (dt_s > 0.1 || dt_s < 0) {
		is_initialized_ = false;
		return;
	}

	// Prediction
	Prediction(dt_s);
	// Update
	if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
	else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar(meas_package);
	}
	// cout << "   Done ProcessMeasurement" << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
	// cout << "Prediction" << endl;
	double dt2 = delta_t * delta_t;
	// cout << "   dt2 = " << dt2 << endl;

	//create example sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented covariance matrix
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	// cout << "      P_aug = " << P_aug << endl;

	//calculate square root of P
	MatrixXd A_aug = P_aug.llt().matrixL();
	// cout << "      A_aug = " << A_aug << endl;

	// Fill x_aug
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;

	//set first column of sigma point matrix
	Xsig_aug.fill(0.0);
	Xsig_aug.col(0) = x_aug;
	//set remaining sigma points
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(i);
	}

	Xsig_pred_.fill(0.0);
	//predict sigma points
	for (int i = 0; i< 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a* dt2 * cos(yaw);
		py_p = py_p + 0.5*nu_a* dt2 * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*dt2;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

	//predicted state mean
	x_.fill(0.0);
	x_ = Xsig_pred_ * weights_;

	//predicted state covariance matrix
	P_.fill(0.0);
	// cout << "      x_ = " << x_.transpose() << endl;
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
											   // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		// // cout << "      Xsig_pred_.col(i) = " << Xsig_pred_.col(i).transpose() << endl;
		//angle normalization
		tools_.NormAngle(x_diff(3));
		tools_.NormAngle(x_diff(4));

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}
	// cout << "      P_ = " << P_ << endl;
	// cout << "   End prediction" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	// cout << "UpdateLidar" << endl;
	int n_z = 2;
	MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);
	Update(meas_package, Zsig);
	// cout << "   End UpdateLidar" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	// cout << "Update Radar" << endl;
	int n_z = 3;
	// predict measurement
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);          //r
		Zsig(1, i) = atan2(p_y, p_x);                  //phi
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}
	Update(meas_package, Zsig);
	// cout << "   Zsig = " << Zsig << endl;
	// cout << "   End Update Radar" << endl;

}

void UKF::Update(MeasurementPackage &meas_package, MatrixXd &Zsig) {
	// cout << "Update()" << endl;
	// cout << "   weights_ = " << weights_.transpose() << endl;
	// cout << "   Zsig = " << Zsig << endl;
	int n_z = Zsig.rows();
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	z_pred = Zsig * weights_;

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											    //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			tools_.NormAngle(z_diff(1));
		}
		// cout << "      z_diff = " << z_diff.transpose() << endl;
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}
	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		R << std_radr_*std_radr_, 0, 0,
			0, std_radphi_*std_radphi_, 0,
			0, 0, std_radrd_*std_radrd_;
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		R << std_laspx_ * std_laspx_, 0,
			0, std_laspy_ * std_laspy_;
	}
	// cout << "   R = " << R << endl;
	S = S + R;

	// cout << "   S = " << S << endl;
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//measurement residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			tools_.NormAngle(z_diff(1));
		}
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		tools_.NormAngle(x_diff(3));

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	// Kalman gain K;
	// Avoid using S.inverse() for better numerical stability
	//K = Tc * S.inverse();
	MatrixXd St = S.transpose();
	MatrixXd Tct = Tc.transpose();
	MatrixXd Kt = St.colPivHouseholderQr().solve(Tct);
	MatrixXd K = Kt.transpose();

	// Residual
	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred;

	//angle normalization
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { // Radar
		tools_.NormAngle(z_diff(1));
	}

	//update state mean and covariance matrix
	// cout << "   K = " << K << endl;
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { // Radar
		NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
		cout << "NIS_radar_ = " << NIS_radar_ << endl;
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) { // Lidar
		NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
		cout << "NIS_lidar_ = " << NIS_laser_ << endl;
	}

	// cout << "   End Update()" << endl;

}
