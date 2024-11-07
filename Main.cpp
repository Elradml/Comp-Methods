/* A program which is intended to examine the application of numerical 
solutions to partial differencial equations using the linear advection problem
* using the explicit and implicit upwind, Lax-Wendroff, and Ritchmyer multi-step methods
using two sets of initial/boundary conditions */
#include <iostream>
#include <vector>
#include <cmath>
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
	return ((analytical_1(x + dx, t) - analytical_1(x, t)) / dt) + (1.75 * (analytical_1(x, t) - analytical_1(x - dx, t)) / dx);
}

//Explicit upwind for set 2
double Eupwind2(double x) {
	return pow(0.5, pow(-x, 2));
}

//Set 1 initial condition
double set1_initial(double x) {
	return (signum(x) + 1) / 2;
}


//Main funciton
int main() {
	//defining vectors
	vector<double> Eupwind_set1;
	vector<double> analytical;
	vector<double> error_v;
	//defining domain
	int min_x = -50;
	int max_x = 50;
	// User enters desired grid spacing
	cout << "Please enter the number of space grid points desired: ";
	double N;
	cin >> N;
	//Set delta t and total time
	double dt = 0.1;
	int max_t = 5;
	//other variables necessary for nested loop
	double solution;
	double x;
	double dx;
	double analytical_solution;
	double error;

	//Nested loops which will calculate the value of the function at each point x
	for (double t = 0; t < max_t - dt; t = t + dt) {
		for (int i = 0; i < N; i++) {
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

			//cout << "x = " << x << "\t\t" << i << "\t\t" << t << endl;
		}
	}
	for (int n = 0; n < Eupwind_set1.size(); n++) {
		cout << "Point " << n << "\t" << Eupwind_set1[n] << "\t" << "Anal\t" << analytical[n] << "\t" << "error = " << error_v[n] << endl;
	}
	return 0;
}