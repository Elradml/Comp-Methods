/* A program which is intended to examine the application of numerical 
solutions to partial differencial equations using the linear advection problem
* using the explicit and implicit upwind, Lax-Wendroff, and Ritchmyer multi-step methods
using two sets of initial/boundary conditions */
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

//Signum function
double signum(double x) {
	if (x > 0) {
		return 1;
	}
	else if (x < 0) {
		return -1;
	}
	else {
		return 0;
	}
}

//Analytical solution for set 1
double analytical_1(double x, double t) {
	return (signum(x - 1.75 * t) + 1) / 2;
}

//Analytical solution for set 2
double analytical_2(double x, double t) {
	return pow(0.5, -pow(x - 1.75 * t, 2));
}

//Explicit upwind for set 1
double Eupwind1(double x, double dx, double dt, double t) {
	return ((analytical_1(x, t + dt) - analytical_1(x, t)) / dt) + (1.75 * (analytical_1(x, t) - analytical_1(x - dx, t)) / dx);
}

//Explicit upwind for set 2
double Eupwind2(double x) {
	return pow(0.5, pow(-x, 2));
}

//Set 1 initial condition
double set1_initial(double x) {
	return (signum(x) + 1) / 2;
}

// Initial conditions for Richtmyer
void initial1(vector<double>& r, double dx) {
	for (int i = 0; i < r.size(); i++) {
		double x = i * dx;
		set1_initial(x);
	}
}

// Boundary conditions for Richtmyer
void boundary1(vector<double>& r, double N) {
	r[0] = 0;
	r[N] = 1;
}

//Main funciton
int main() {
	// ------------------------
	// Explicit Upwind Solution
	// ------------------------
	
	//defining vectors
	vector<double> Eupwind_set1;
	vector<double> analytical;
	vector<double> error_v;
	//defining domain
	const int min_x = -50;
	const int max_x = 50;
	// User enters desired grid spacing
	cout << "Please enter the number of space grid points desired: ";
	double N;
	cin >> N;
	//Set delta t and total time
	double dt = 0.1;
	int max_t = 2;
	//other variables necessary for nested loop
	double solution;
	double x;
	double dx;
	double analytical_solution;
	double error;

	//Time step loop
	for (double t = 0; t < max_t - dt; t = t + dt) {
		// Spacial step loop
		for (double i = 0; i < N; i++) {
			x = min_x + (max_x - min_x) * (i / (N - 1));
			analytical_solution = analytical_1(x, t);
			if (i == 1) {
				dx = x - min_x;
			}
			if (t == 0) {
				solution = set1_initial(x);
			}
			else {
				solution = Eupwind1(x, dx, dt, t);
			}

			//Placing solutions into vectors
			Eupwind_set1.push_back(solution);
			analytical.push_back(analytical_solution);
			error = solution - analytical_solution;
			error_v.push_back(error);
		}
	}
	/*for (int n = 0; n < Eupwind_set1.size(); n++) {
		cout << "Point " << n << "\t" << Eupwind_set1[n] << "\t" << "Anal\t" << analytical[n] << "\t" << "error = " << error_v[n] << endl;
	}*/

	// ------------------------
	// Implicit Upwind Solution
	// ------------------------
	
	// Important Variables
	const double u = 1.75;
	double courant = (u * dt) / dx;

	// Initialize Thomas algorithm vectors & vectors for solutions
	vector<double> f(N, 0.0);
	vector<double> fn(N, 0.0);
	vector<double> sub_diag(N, -courant);
	vector<double> main_diag(N, 1 + courant);
	vector<double> sup_diag(N, 0.0);

	// Set initial conditions
	for (int i = 0; i < N-1; i++) {
		f[i] = set1_initial(i);
	}

	// Time step loop
	for (double n = 0; n < max_t; n = n + dt) {
		//vector<double> super(N, 0.0);
		vector<double> fin(N, 0.0);
		
		// Boundary conditions
		fin[0] = 0;
		fin[N-1] = 1;

		//cout << "f: " << f.size() << endl << "fn: " << fn.size() << endl << "sub: " << sub_diag.size() << endl << "main: " << main_diag.size() << endl << "sup: " << sup_diag.size() << endl << "fin: " << fin.size() << endl;

		// Forward Sweep
		for (double i = 1; i < N; i++) {
			double m = sub_diag[i] / main_diag[i - 1];
			main_diag[i] = main_diag[i] - m * sup_diag[i - 1];
			fin[i] = f[i] - m * fin[i - 1];
		}

		// Backward Substitution
		fn[N - 1] = fin[N - 1] / main_diag[N - 1];
		for (double i = N - 2; i >= 0; i--) {
			fn[i] = (fin[i] - sup_diag[i] * fn[i + 1]) / main_diag[i];
		}

		f = fn;

		/*if (fmod(n, 10) == 0) {
			cout << "Time Step " << n << endl;
			for (const auto& val : f) {
				cout << val << " ";
			}
			cout << endl;
		}
		for (int i = 0; i < f.size(); i++) {
			cout << i << "\t" << f[i] << endl;
		}*/
	}

	// -----------------------------
	// Richtmyer Multi-Step Solution
	// -----------------------------
	
	// Vectors for solutions
	vector<double> r(N + 2, 0);
	vector<double> r_half(N + 2, 0);
	vector<double> r_new(N + 2, 0);

	initial1(r, dx);

	// Time step loop
	for (double n = 0; n < N; n++) {
		// Apply boundaries
		boundary1(r, N);

		// Predictor step
		for (int i = 1; i < N + 1; i++) {
			r_half[i] = 0.5 * (r[i] + r[i + 1]) - 0.5 * courant * (r[i + 1] - r[i]);
		}

		r_half[0] = 0;
		r_half[N + 1] = 1;

		// Corrector step
		for (int i = 1; i < N; i++) {
			r_new[i] = r[i] - courant * (r_half[i] - r_half[i - 1]);
		}
		r = r_new;

		if (fmod(n, 10) == 0) {
			cout << "Time step " << n << ":" << endl;
			for (int i = 1; i < N + 1; i++) {
				cout << r[i] << " ";
			}
			cout << endl;
		}
	}
	
	return 0;
}