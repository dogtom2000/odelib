#include "Integrate.h"
#include <algorithm>
#include <cmath>

// logic to use switch statement as an arithmetic if
#define aif(a) ((a > 0) ? 'p' : ((a < 0) ? 'n' : '0'))
// plus 1, minus 1 used to convert between 0/1 based indices
#define p1(a) (a + 1)
#define m1(a) (a - 1)
// 2d addressing logic for phi
#define phi(j, i) phi[j + i * neqn]

Integrate::Integrate()
{
}

Integrate::~Integrate()
{
}

void Integrate::take_step()
{
	// test input parameters, initialize
	// if input is invalid return
	block0();
	if (crash) { return; }

	while (true)
	{
		// compute coefficients of formulas
		block1();
		// predict solution, estimate local errors
		block2();

		if (step_fail)
		{
			// failed step:
			// restore input parameters, select new order/step size
			// if step size is too small return, otherwise go to block 1
			block3();
			if (crash) { return; }
		}
		else
		{
			// successful step:
			// correct solution, select new order/step size, return
			block4();
			return;
		}
	}
}

void Integrate::interp()
{
	double hi = xout - x;
	unsigned int ki = kold + 1;
	unsigned int kip1 = ki + 1;

	for (size_t i = 0; i < ki; i++)
	{
		size_t temp1 = i + 1;
		wi[i] = 1.0 / temp1;
	}

	double term = 0.0;

	for (size_t j = 1; j < ki; j++)
	{
		size_t jm1 = j - 1;
		double psijm1 = psi[jm1];
		double gamma = (hi + term) / psijm1;
		double eta = hi / psijm1;
		size_t limit1 = kip1 - j - 1;
		for (size_t i = 0; i < limit1; i++)
		{
			wi[i] = gamma * wi[i] - eta * wi[i + 1];
		}
		gi[j] = wi[0];
		rho[j] = gamma * rho[jm1];
		term = psijm1;
	}

	for (size_t l = 0; l < neqn; l++)
	{
		yout[l] = 0.0;
		ypout[l] = 0.0;
	}

	for (size_t j = 0; j < ki; j++)
	{
		size_t i = kip1 - p1(j);
		double temp2 = gi[m1(i)];
		double temp3 = rho[m1(i)];
		for (size_t l = 0; l < neqn; l++)
		{
			yout[l] += temp2 * phi(l, m1(i));
			ypout[l] += temp3 * phi(l, m1(i));
		}
	}
	for (size_t l = 0; l < neqn; l++)
	{
		yout[l] = y[l] + hi * yout[l];
	}
}

void Integrate::extrap()
{
	f(x, y, yp);
	for (size_t l = 0; l < neqn; l++)
	{
		yout[l] = y[l] + h * yp[l];
	}
}

void Integrate::block0()
{
	// test input parameters, initialize
	// if step size or epsilon are too small crash and return
	crash = false;
	test_inputs();
	if (crash) { return; }

	// if this is the first step inialize starting parameters
	if (start) { initialize(); }

	// initialize failed step counter to 0
	ifail = 0;
}

void Integrate::block1()
{
	// set offset values for k
	kp1 = k + 1;
	kp2 = k + 2;
	km1 = k - 1;
	km2 = k - 2;

	// if step size has changed reset step size counter
	if (h != hold) { ns = 0; }

	// increase step size counter, limit to kold + 1
	ns = std::min(ns + 1, kold + 1);

	// coefficients are only updated if k >= ns
	if (k >= ns)
	{
		compute_coefficients();

		if (ns == 1) { initialize_vw(); }
		else { update_vw(); }

		if (k > ns) { compute_g(); }
	}
}

void Integrate::block2()
{
	// phistar is only calculated if k > ns
	if (k > ns) { phi_star(); }

	// predict the value p at x + h
	predict();

	// increment x and calculate derivatives at x + h, p
	xold = x;
	x += h;
	f(x, p, yp);

	// estimate the error at orders k, k-1, k-2
	// if error is greater than epsilon the step fails
	estimate_error();
	step_fail = err <= eps ? false : true;
}

void Integrate::block3()
{
	// on a step failure end the starting phase
	phase1 = false;

	// restore input parameters
	restore();

	// increment failed step counter
	ifail++;

	// halve step size
	// if more than 2 failed steps reset step size
	// if more than 3 failed steps reduce order to 1
	update_o_h_fail();
}

void Integrate::block4()
{
	kold = k;
	hold = h;

	// correct the value p to y at n + 1
	correct();

	// calculate the derivatives at x(n + 1), y(n + 1)
	f(x, y, yp);

	// update the divided differences (phi) for the next step
	update_differences();

	// if k was reduced or k = 12 the end the starting phase
	// if in the starting phase increase the order and double step size
	// otherwise determine new order and step size
	erkp1 = 0.0;
	if ((knew == km1) || (k == 12)) { phase1 = false; }
	if (phase1)
	{
		k = kp1;
		erk = erkp1;
		h *= 2;
	}
	else
	{
		update_o_h_success();
	}
}

void Integrate::test_inputs()
{
	// if step size is too small, increase step size, crash, and return
	if (std::abs(h) < fouru * std::abs(x))
	{
		h = fouru * std::abs(x) * sgn(h);
		crash = true;
		return;
	}

	// calculate round, round is greater than the largest value in y multiplied by 2u
	round = 0.0;
	for (size_t l = 0; l < neqn; l++)
	{
		round += pow(y[l] / wt[l], 2);
	}
	round = twou * sqrt(round);

	// if epsilon is too small, increase epsilon, crash, and return
	if (0.5 * eps < round)
	{
		eps = 2.0 * round * (1.0 + fouru);
		crash = true;
		return;
	}
}

void Integrate::initialize()
{
	// initialize values
	hold = 0.0;
	k = 1;
	kold = 0;
	start = false;
	phase1 = true;
	nornd = true;

	// evaluate derivatives at x, y
	f(x, y, yp);

	// initialize phi
	// calculate sum, sum is greater than the largeset value of yp
	double sum = 0.0;
	for (size_t l = 0; l < neqn; l++)
	{
		phi(l, 0) = yp[l];
		phi(l, 1) = 0.0;
		sum += pow(yp[l] / wt[l], 2);
	}
	sum = sqrt(sum);

	// if step size is too small for desired epsilon, increase step size
	absh = std::abs(h);
	if (eps < 16.0 * sum * h * h)
	{
		absh = 0.25 * sqrt(eps / sum);
	}
	h = std::max(absh, fouru * std::abs(x)) * sgn(h);

	// if epsilon is less than 100 * round, use propogated round off control
	if (0.5 * eps < 100.0 * round)
	{
		nornd = false;
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, 14) = 0.0;
		}
	}
}

void Integrate::compute_coefficients()
{
	beta[m1(ns)] = 1.0;
	alpha[m1(ns)] = 1.0 / ns;
	double temp1 = h * ns;
	sigma[ns] = 1.0;
	if (k > ns)
	{
		for (size_t i = ns; i < k; i++)
		{
			size_t im1 = i - 1;
			double temp2 = psi[im1];
			psi[im1] = temp1;
			beta[i] = beta[im1] * psi[im1] / temp2;
			temp1 = temp2 + h;
			alpha[i] = h / temp1;
			sigma[i + 1] = p1(i) * alpha[i] * sigma[i];
		}
	}
	psi[m1(k)] = temp1;
}

void Integrate::initialize_vw()
{
	for (size_t iq = 0; iq < k; iq++)
	{
		double temp3 = p1(iq) * (p1(iq) + 1);
		v[iq] = 1.0 / temp3;
		w[iq] = v[iq];
	}
}

void Integrate::update_vw()
{
	if (k > kold)
	{
		double temp4 = k * kp1;
		v[m1(k)] = 1.0 / temp4;
		if (ns >= 3)
		{
			for (size_t j = 0; j < ns - 2; j++)
			{
				size_t i = k - p1(j);
				v[m1(i)] -= alpha[j + 1] * v[m1(i) + 1];
			}
		}
	}

	size_t limit1 = kp1 - ns;
	double temp5 = alpha[m1(ns)];
	for (size_t iq = 0; iq < limit1; iq++)
	{
		v[iq] -= temp5 * v[iq + 1];
		w[iq] = v[iq];
	}
	g[ns] = w[0];
}

void Integrate::compute_g()
{
	for (size_t i = ns + 1; i < kp1; i++)
	{
		size_t limit2 = kp2 - p1(i);
		double temp6 = alpha[i - 1];
		for (size_t iq = 0; iq < limit2; iq++)
		{
			w[iq] -= temp6 * w[iq + 1];
		}
		g[i] = w[0];
	}
}

void Integrate::phi_star()
{
	for (size_t i = ns; i < k; i++)
	{
		double temp1 = beta[i];
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, i) *= temp1;
		}
	}
}

void Integrate::predict()
{
	for (size_t l = 0; l < neqn; l++)
	{
		phi(l, m1(kp2)) = phi(l, m1(kp1));
		phi(l, m1(kp1)) = 0.0;
		p[l] = 0.0;
	}

	for (size_t j = 0; j < k; j++)
	{
		size_t i = m1(kp1 - p1(j));
		size_t ip1 = i + 1;
		double temp2 = g[i];
		for (size_t l = 0; l < neqn; l++)
		{
			p[l] += temp2 * phi(l, i);
			phi(l, i) += phi(l, ip1);
		}
	}

	if (nornd)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			p[l] = y[l] + h * p[l];
		}
	}
	else
	{
		for (size_t l = 0; l < neqn; l++)
		{
			double tau = h * p[l] - phi(l, 14);
			p[l] = y[l] + tau;
			phi(l, 15) = (p[l] - y[l]) - tau;
		}
	}
}

void Integrate::estimate_error()
{
	absh = std::abs(h);
	erkm2 = 0.0;
	erkm1 = 0.0;
	erk = 0.0;

	for (size_t l = 0; l < neqn; l++)
	{
		double temp3 = 1.0 / wt[l];
		double temp4 = yp[l] - phi(l, 0);

		switch (aif(km2))
		{
		case 'p':
			erkm2 += pow((phi(l, m1(km1)) + temp4) * temp3, 2);
		case '0':
			erkm1 += pow((phi(l, m1(k)) + temp4) * temp3, 2);
		case 'n':
			erk += pow(temp4 * temp3, 2);
		}
	}

	switch (aif(km2))
	{
	case 'p':
		erkm2 = absh * sigma[m1(km1)] * gstar[m1(km2)] * sqrt(erkm2);
	case '0':
		erkm1 = absh * sigma[m1(k)] * gstar[m1(km1)] * sqrt(erkm1);
	case 'n':
		break;
	}

	double temp5 = absh * sqrt(erk);
	double err = temp5 * (g[m1(k)] - g[m1(kp1)]);
	erk = temp5 * sigma[m1(kp1)] * gstar[m1(k)];
	knew = k;

	switch (aif(km2))
	{
	case 'p':
		if (std::max(erkm1, erkm2) <= erk)
		{
			knew = km1;
		}
		break;
	case '0':
		if (erkm1 <= 0.5 * erk)
		{
			knew = km1;
		}
		break;
	case 'n':
		break;
	}
}

void Integrate::restore()
{
	x = xold;

	for (size_t i = 0; i < k; i++)
	{
		double temp1 = 1.0 / beta[i];
		size_t ip1 = i + 1;
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, i) = temp1 * (phi(l, i) - phi(l, ip1));
		}
		if (k > 1)
		{
			for (size_t i = 1; i < k; i++)
			{
				psi[i - 1] = psi[i] - h;
			}
		}
	}
}

void Integrate::update_o_h_fail()
{
	double temp2 = 0.5;
	if (ifail > 3)
	{
		if (0.5 * eps < 0.25 * erk)
		{
			temp2 = sqrt(0.5 * eps / erk);
		}
	}
	if (ifail >= 3)
	{
		knew = 1;
	}
	h *= temp2;
	k = knew;
	if (std::abs(h) < fouru * std::abs(x))
	{
		crash = true;
		h = fouru * std::abs(x) * sgn(h);
		eps *= 2;
	}
}

void Integrate::correct()
{
	double temp1 = h * g[m1(kp1)];
	if (nornd)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			y[l] = p[l] + temp1 * (yp[l] - phi(l, 0));
		}
	}
	else
	{
		for (size_t l = 0; l < neqn; l++)
		{
			double rho = temp1 * (yp[l] - phi(l, 0)) - phi(l, 15);
			y[l] = p[l] + rho;
			phi(l, 14) = (y[l] - p[l]) - rho;
		}
	}
}

void Integrate::update_differences()
{
	for (size_t l = 0; l < neqn; l++)
	{
		phi(l, m1(kp1)) = yp[l] - phi(l, 0);
		phi(l, m1(kp2)) = phi(l, m1(kp1)) - phi(l, m1(kp2));
	}
	for (size_t i = 0; i < k; i++)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, i) += phi(l, m1(kp1));
		}
	}
}

void Integrate::update_o_h_success()
{
	if (knew == km1)
	{
		k = km1;
		erk = erkm1;
	}
	else if (k < ns)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			erkp1 += pow(phi(l, m1(kp2)) / wt[l], 2);
		}
		erkp1 = absh * gstar[m1(kp1)] * sqrt(erkp1);

		if (k == 1)
		{
			if (erkp1 < 0.5 * erk)
			{
				k = kp1;
				erk = erkp1;
			}
		}
		else
		{
			if (erkm1 <= std::min(erk, erkp1))
			{
				k = km1;
				erk = erkm1;
			}
			else if ((erkp1 < erk) && (k < 12))
			{
				k = kp1;
				erk = erkp1;
			}
		}
	}

	double hnew = h + h;
	if (0.5 * eps >= erk * two[k])
	{
		h = hnew;
		return;
	}
	hnew = h;
	if (0.5 * eps >= erk)
	{
		h = hnew;
		return;
	}
	size_t temp2 = k + 1;
	double r = pow((0.5 * eps / erk), 1.0 / temp2);
	hnew = absh * std::max(0.5, std::min(0.9, r));
	hnew = std::max(hnew, fouru * std::abs(x)) * sgn(h);
	h = hnew;
	return;
}