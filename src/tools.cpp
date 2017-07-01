#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	for(int i=0; i < estimations.size(); ++i)
	{
		VectorXd temp = (estimations[i] - ground_truth[i]);
		temp = temp.array()*temp.array();
		rmse += temp;

	}

	rmse = rmse/estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);


	float mod = px*px + py*py;
	float rootMod = sqrt(mod);

	//check division by zero
	if(fabs(mod) < 0.0001)
	{
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	Hj << px/rootMod, py/rootMod, 0, 0,
			-py/mod, px/mod, 0, 0,
			(py*(vx*py-vy*px))/sqrt(mod*mod*mod),(px*(vx*py-vy*px))/sqrt(mod*mod*mod), px/rootMod, py/rootMod;

	return Hj;
}
