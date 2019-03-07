#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
  			  0, 1, 0, 0;
  
   // Initial transition matrix F_
      ekf_.F_ = MatrixXd(4, 4);
      ekf_.F_ << 1, 0, 1, 0,
                 0, 1, 0, 1,
                 0, 0, 1, 0,
                 0, 0, 0, 1;
    
    //Initial state covariance matrix P_
      ekf_.P_ = MatrixXd(4, 4);
      ekf_.P_ << 1, 0, 0, 0,
            	 0, 1, 0, 0,
                 0, 0, 1000, 0,
                 0, 0, 0, 1000;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {

    // Initialize the state ekf_.x_ with the first measurement.
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      //First measurement of RADAR sensor
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rhodot = measurement_pack.raw_measurements_[2];
      
	  // Conversion of radar measurements from polar to cartesian coordinates 
      double px = phi * sin(rho);
      double py = phi * cos(rho);
      double vx = rhodot * sin(rho);
      double vy = rhodot * cos(rho);
      
      //Initial state for Radar Measurements
	  ekf_.x_ << px, py, vx, vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    //First measurement of LASER sensor
      
      //Initialize state for Laser Measurements.
      double px =  measurement_pack.raw_measurements_[0];
      double py =  measurement_pack.raw_measurements_[1];
      double vx =  0;
      double vy =  0;
      
      ekf_.x_ << px, py, vx, vy;
    }

    
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }
   /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //Calculating the elapsed time between current and previous timestamp
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  //Transition F matrix integrated with time
  ekf_.F_(0, 2) = dt;   //0 row and 2 column
  ekf_.F_(1, 3) = dt;   //1 row and 3 column
  
  //set up noise component
  float noise_ax = 9;
  float noise_ay = 9;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  //Updation of process covariance matrix Q according to new elapsed time.
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  //Prediction of position and velocity of bicycle
  ekf_.Predict();

  /*****************************************************************************
   *  Updation
   ****************************************************************************/
  //Updation of position and velocity of bicycle according to new measuremnets for Laser and Radar
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
  cout << "EKF : First measurement RADAR" << endl;
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);    
  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the predicted state x_ and covariance Matrix P_
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
