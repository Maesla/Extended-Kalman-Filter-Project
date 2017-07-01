#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.14159265

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::UpdateF(float dt)
{
    	//state covariance matrix P
	F_ = MatrixXd(4, 4);
	F_ << 1, 0, dt, 0,
		  0, 1, 0, dt,
		  0, 0, 1, 0,
		  0, 0, 0, 1;
}

void KalmanFilter::UpdateQ(float dt, float noise_ax, float noise_ay)
{
    float dt2 = dt*dt;
    float dt3 = dt2*dt;
    float dt4 = dt3*dt;

    float c00 = 0.25f*dt4*noise_ax;
    float c02 = 0.5f*dt3*noise_ax;

    float c11 = 0.25f*dt4*noise_ay;
    float c13 = 0.5*dt3*noise_ay;

    float c20 = 0.5*dt3*noise_ax;
    float c22 = dt2*noise_ax;

    float c31 = 0.5*dt3*noise_ay;
    float c33 = dt2*noise_ay;

    Q_ = MatrixXd(4, 4);
	Q_ << c00, 0, c02, 0,
		  0, c11, 0, c13,
		  c20, 0, c22, 0,
		  0, c31, 0, c33;
}

void KalmanFilter::Predict()
{
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	CalculateJacobian();
	VectorXd z_pred = TransformStateFromCartesian2Polar();
	VectorXd y = z - z_pred;
	y = KeepRhoInRangeMinusPiAndPi(y);
	MatrixXd Ht = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;
}

void KalmanFilter::CalculateJacobian()
{
	Hj_ = MatrixXd(3, 4);
	//recover state parameters
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);


	float mod = px*px + py*py;
	float rootMod = sqrt(mod);

	//check division by zero
	if(fabs(mod) < 0.0001){
		//cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		//return Hj;
		return;
	}

	Hj_ << px/rootMod, py/rootMod, 0, 0,
		-py/mod, px/mod, 0, 0,
		(py*(vx*py-vy*px))/sqrt(mod*mod*mod),(px*(vx*py-vy*px))/sqrt(mod*mod*mod), px/rootMod, py/rootMod;
	//TODO: YOUR CODE HERE

	//check division by zero

	//compute the Jacobian matrix
}

Eigen::VectorXd KalmanFilter::TransformStateFromCartesian2Polar()
{
	VectorXd polar(3);
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	float polar_0 = sqrt(px*px + py*py);
	float polar_1 = atan2(py,px);
	float polar_2 = (px*vx + py*vy)/polar_0;

	polar << polar_0, polar_1, polar_2;

	return polar;
}

Eigen::VectorXd KalmanFilter::KeepRhoInRangeMinusPiAndPi(const VectorXd &y)
{
	VectorXd new_y(3);
	float rho = y(1);

	while (rho < -PI || rho > PI)
	{
		if(rho < -PI)
		{
			rho = rho + 2*PI;
		}
		else if (rho > PI)
		{
			rho = rho - 2*PI;
		}
	}

	new_y << y(0), rho, y(2);
	return new_y;
}


