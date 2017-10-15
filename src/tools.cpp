#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0.0, 0.0, 0.0, 0.0;

	if (estimations.size() == 0) {
		cout << "Estimation is emptys" << endl;
		return rmse;
	}

	// The estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size()) {
		cout << "Size of estimations and ground_truth must match" << endl;
		return rmse;
	}

	// Accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		VectorXd residual = estimations[i] - ground_truth[i];

		// Coefficient-wise multiplication
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// Calculate the mean
	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();
	// cout << "RMSE = " << rmse.transpose() << endl << endl;

	return rmse;
}

void Tools::NormAngle(double& input) {
	//while (input > M_PI) input -= 2.*M_PI;
	//while (input < -M_PI) input += 2.*M_PI;
	input = atan2(sin(input), cos(input));
	return;
}

