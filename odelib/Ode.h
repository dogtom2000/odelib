#pragma once

#include "Integrate.h"

class Ode
{
public:
	// constructor and destructor
	Ode(void(*)(double, double*, double*), unsigned int, double*, double, double, double, double, double*);
	~Ode();

	// logical
	bool stiff;

	// variables
	int iflag;
	int isgn;
	int isgnold;
	int delsgn;

	unsigned int neqn;
	unsigned int nostep;
	unsigned int maxnum{ 5000 };
	unsigned int kle4;

	double relerr;
	double abserr;
	double releps;
	double abseps;

	double t;
	double tout;
	double tend;

	double del;
	double absdel;

	// arrays
	double* y;

	// integrator
	Integrate integrate;

	// member functions
	void step();
	bool test_inputs();
	void setup();
	void first_step();
	void set_weights();
	void increment();
	void end_interp();
	void end_extrap();
	void end_work();
	void end_tol();

	// inline functions
	template <class T>
	int sgn(T a) { return (a > T(0)) - (a < T(0)); }
};