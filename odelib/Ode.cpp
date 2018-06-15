#include "Ode.h"
#include <algorithm>
#include <cmath>

Ode::Ode(void(*f)(double, double*, double*), unsigned int neqn, double* y, double t, double tout, double relerr, double abserr, double* work)
	: neqn(neqn), y(y), t(t), tout(tout), relerr(relerr), abserr(abserr)
{
	iflag = 1;

	integrate.y = work + 0 * neqn;
	integrate.yp = work + 1 * neqn;
	integrate.yout = work + 2 * neqn;
	integrate.ypout = work + 3 * neqn;
	integrate.p = work + 4 * neqn;
	integrate.wt = work + 5 * neqn;
	integrate.phi = work + 6 * neqn;

	integrate.neqn = neqn;
	integrate.f = f;

	double halfu = 0.5;
	while (1 + halfu > 1)
	{
		halfu *= 0.5;
	}
	integrate.twou = 4 * halfu;
	integrate.fouru = 8 * halfu;
}

Ode::~Ode()
{
}

void Ode::step()
{
	if (!test_inputs())
	{
		return;
	}

	setup();

	if ((iflag == 1) || (isgnold < 0) || (delsgn * del <= 0))
	{
		first_step();
	}

	while (nostep < maxnum)
	{
		if (std::abs(integrate.x - t) >= absdel)
		{
			end_interp();
			return;
		}

		if ((isgn < 0) && (std::abs(tout - integrate.x) < integrate.fouru * std::abs(integrate.x)))
		{
			end_extrap();
			return;
		}

		set_weights();

		integrate.take_step();

		if (integrate.crash)
		{
			end_tol();
			return;
		}
		else
		{
			increment();
		}
	}

	end_work();
	return;
}

bool Ode::test_inputs()
{
	if (t == tout) { iflag = 6; return; }
	if (relerr < 0 || abserr < 0) { iflag = 6; return; }
	integrate.eps = std::max(relerr, abserr);
	isgn = sgn(iflag);
	iflag = std::abs(iflag);
	return true;
}

void Ode::setup()
{
	del = tout - t;
	absdel = std::abs(del);
	tend = t + 10.0 * del;
	if (isgn < 0) { tend = tout; }
	nostep = 0;
	kle4 = 0;
	stiff = false;
	releps = relerr / integrate.eps;
	abseps = abserr / integrate.eps;
}

void Ode::first_step()
{
	integrate.start = true;
	integrate.x = t;
	delsgn = sgn(del);
	integrate.h = sgn(tout - integrate.x) * std::max(std::abs(tout - integrate.x), integrate.fouru * std::abs(integrate.x));
	for (size_t l = 0; l < neqn; l++)
	{
		integrate.y[l] = y[l];
	}
}

void Ode::set_weights()
{
	integrate.h = sgn(integrate.h) * std::min(std::abs(integrate.h), std::abs(tend - integrate.x));
	for (size_t l = 0; l < neqn; l++)
	{
		integrate.wt[l] = releps * std::abs(integrate.y[l]) + abseps;
	}
}

void Ode::increment()
{
	nostep++;
	kle4++;
	if (integrate.kold > 4) { kle4 = 0; }
	if (kle4 > 50) { stiff = true; }
}

void Ode::end_interp()
{
	integrate.xout = tout;
	integrate.interp();
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = integrate.yout[l];
	}
	iflag = 2;
	t = tout;
	isgnold = isgn;
}

void Ode::end_extrap()
{
	integrate.h = tout - integrate.x;
	integrate.extrap();
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = integrate.yout[l];
	}
	iflag = 2;
	t = tout;
	isgnold = isgn;
}

void Ode::end_work()
{
	iflag = isgn * 4;
	if (stiff) { iflag = isgn * 5; }
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = integrate.y[l];
	}
	t = integrate.x;
	isgnold = 1;
}

void Ode::end_tol()
{
	iflag = isgn * 3;
	relerr = integrate.eps * releps;
	abserr = integrate.eps * abseps;
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = integrate.y[l];
	}
	t = integrate.x;
	isgnold = 1;
}