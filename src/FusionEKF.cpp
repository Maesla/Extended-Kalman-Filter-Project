#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
  			  0, 1, 0, 0,
  			  0, 0, 1, 0,
  			  0, 0, 0, 1;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
  			  0, 1, 0, 0;

  noise_ax = 9.0;
  noise_ay = 9.0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;

	ekf_.P_ = P_;
	ekf_.H_ = H_laser_;
	 previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    	//cout << "RADAR" << endl;
    	/**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

    	ekf_.x_ = VectorXd(4);
    	float phi = measurement_pack.raw_measurements_[0];
    	float rho = measurement_pack.raw_measurements_[1];
    	float px = cos(rho)*phi;
    	float py = sin(rho)*phi;
    	ekf_.x_  << px, py, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
    	//cout << "LASER" << endl;

    	ekf_.x_ = VectorXd(4);
    	ekf_.x_  << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;


      /**
      Initialize state.
      */
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  	previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  ekf_.UpdateF(dt);
  ekf_.UpdateQ(dt, noise_ax, noise_ay);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
  	//cout << "RADAR" << endl;
	ekf_.R_ = R_radar_;

	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    // Radar updates
  } else {
    // Laser updates
	 // cout << "LASER" << endl;
	ekf_.R_ = R_laser_;

	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
