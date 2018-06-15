#include "stdafx.h"
#include <cmath>

namespace diffeq {

	void two_body(double x, double* y, double* yp)
	{
		double temp = y[0] * y[0] + y[1] * y[1];
		double mu = 1;
		temp *= std::sqrt(temp);
		yp[0] = y[2];
		yp[1] = y[3];
		yp[2] = -mu / temp * y[0];
		yp[3] = -mu / temp * y[1];
	}

	void seven_body(double x, double* y, double* yp)
	{
		unsigned int nbody = 7;
		double mu[7] = { 1.0, 0.00000244784, 0.000003003489663, 0.000000036958626, 0.0000003227155, 0.000954594263, 0.00028581485 };
		for (size_t i = 0; i < nbody; i++)
		{
			int ii = 6 * i;
			for (size_t j = 0; j < 3; j++)
			{
				yp[ii + j] = y[ii + j + 3];
			}

			yp[ii + 3] = 0.0;
			yp[ii + 4] = 0.0;
			yp[ii + 5] = 0.0;

			for (size_t j = 0; j < nbody; j++)
			{
				if (i != j)
				{
					int jj = 6 * j;
					double rx = y[ii + 0] - y[jj + 0];
					double ry = y[ii + 1] - y[jj + 1];
					double rz = y[ii + 2] - y[jj + 2];
					double r = std::sqrt(rx * rx + ry * ry + rz * rz);

					yp[ii + 3] -= mu[j] * rx / std::pow(r, 3);
					yp[ii + 4] -= mu[j] * ry / std::pow(r, 3);
					yp[ii + 5] -= mu[j] * rz / std::pow(r, 3);
				}
			}
			if (i != 0)
			{
				yp[ii + 3] -= yp[3];
				yp[ii + 4] -= yp[4];
				yp[ii + 5] -= yp[5];
			}
		}
		yp[3] = 0.0;
		yp[4] = 0.0;
		yp[5] = 0.0;
	}
}