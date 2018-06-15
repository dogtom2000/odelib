#pragma once

class Integrate
{
public:
	// constants
	const double gstar[13]{ 0.500, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679, 0.00592, 0.00524, 0.00468 };
	const double two[13]{ 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192 };

	// constructor and destructor
	Integrate();
	~Integrate();

	// logical
	bool start;
	bool crash;
	bool phase1;
	bool nornd;
	bool step_fail;

	// variables
	unsigned int neqn;
	unsigned int ifail;

	unsigned int ns;
	unsigned int k;
	unsigned int knew;
	unsigned int kold;
	unsigned int kp1;
	unsigned int kp2;
	int km1;
	int km2;

	double twou;
	double fouru;

	double eps;
	double err;
	double round;
	double erk;
	double erkp1;
	double erkm1;
	double erkm2;

	double h;
	double hold;
	double absh;

	double x;
	double xold;
	double xout;

	// variable sized arrays;
	double* y;
	double* yp;
	double* yout;
	double* ypout;
	double* wt;
	double* p;
	double* phi;

	// fixed size arrays
	double alpha[12];
	double beta[12];
	double psi[12];
	double sigma[13]{ 1.0 };

	double v[12];
	double w[12];
	double g[13]{ 1.0, 0.5 };

	double wi[13];
	double gi[13]{ 1.0 };
	double rho[13]{ 1.0 };

	// member functions
	void(*f)(double, double*, double*);

	// functions called from Ode
	void take_step();
	void interp();
	void extrap();

	// block functions
	void block0();
	void block1();
	void block2();
	void block3();
	void block4();

	// block 0 functions
	void test_inputs();
	void initialize();

	// block 1 functions
	void compute_coefficients();
	void initialize_vw();
	void update_vw();
	void compute_g();

	// block 2 functions
	void phi_star();
	void predict();
	void estimate_error();

	// block 3 functions
	void restore();
	void update_o_h_fail();

	// block 4 functions
	void correct();
	void update_differences();
	void update_o_h_success();

	// inline functions
	template <class T>
	int sgn(T a) { return (a > T(0)) - (a < T(0)); }
};