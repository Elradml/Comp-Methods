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
	return pow(0.5, pow(-(x - 1.75 * t), 2));
}

//Explicit upwind for set 1
double Eupwind1(double x, double dx, double dt, double t) {
	return ((analytical_1(x, t + dt) - analytical_1(x, t)) / dt) + 1.75 * (analytical_1(x, t) - analytical_1(x - dx, t) / dx);
}

//Explicit upwind for set 2
double Eupwind2(double x, double dx, double dt, double t) {
	return (analytical_2(x, t + dt) - analytical_2(x, t) / dt) + 1.75 * (analytical_2(x, t) - analytical_2(x - dx, t) / dx);
}

//Set 1 initial condition
double set1_initial(double x) {
	return (signum(x) + 1) / 2;
}

// Set 2 initial condition
double set2_initial(double x) {
	return pow(0.5, pow(-x, 2));
}

// Initial conditions for Richtmyer Set 1
void initial1(vector<double>& r, double dx) {
	for (int i = 0; i < r.size(); i++) {
		double x = i * dx;
		r[i] = set1_initial(x);
	}
}

// Initial conditions for Richtmyer Set 2
void initial2(vector<double>& r, double dx) {
	for (int i = 0; i < r.size(); i++) {
		double x = i * dx;
		r[i] = set2_initial(x);
	}
}

// Boundary conditions for Richtmyer Set 1
void boundary1(vector<double>& r, double N) {
	r[0] = 0;
	r[N] = 1;
}

// Boundary conditions for Richtmyer Set 2
void boundary2(vector<double>& r, double N) {
	r[0] = 0;
	r[N] = 0;
}

//Main funciton
int main() {
	// ------------------------
	// Explicit Upwind Solution
	// ------------------------
	
	//defining vectors
	vector<double> Eupwind_set1;
	vector<double> analytical;
	vector<double> analytical2;
	vector<double> error_v;
	vector<double> x_values;
	//defining domain
	const int min_x = -50;
	const int max_x = 50;
	// User enters desired grid spacing
	cout << "Please enter the number of space grid points desired: ";
	double N;
	cin >> N;
	cout << N << endl;
	//Set delta t and total time
	double dt = 0.1;
	int max_t = 2;
	//other variables necessary for nested loop
	double solution;
	double solution2;
	double x;
	double dx;
	double analytical_solution1;
	double analytical_solution2;
	double error1;
	double error2;
	cout << "***********************" << endl << "Explicit Solutions" << endl << "***********************" << endl;

	//Time step loop
	for (double t = 0; t < max_t + dt; t += dt) {
		// Spacial step loop
		for (double i = 0; i < N; i++) {
			x = min_x + (max_x - min_x) * (i / (N - 1));
			x_values.push_back(x);
			analytical_solution1 = analytical_1(x, t);
			analytical_solution2 = analytical_2(x, t);
			if (i == 1) {
				dx = x - min_x;
			}
			if (t == 0) {
				solution = set1_initial(x);
				solution2 = set2_initial(x);
			}
			else {
				solution = Eupwind1(x, dx, dt, t);
				solution2 = Eupwind2(x, dx, dt, t);
			}

			//Placing solutions into vectors
			/*Eupwind_set1.push_back(solution);
			analytical.push_back(analytical_solution1);
			analytical2.push_back(analytical_solution2);
			error = solution - analytical_solution2;
			error_v.push_back(error);*/
			error1 = solution - analytical_solution1;
			error2 = solution2 - analytical_solution2;

			cout << "t = " << t << "\t" << "x = " << x << "\t" << "A1 = " << analytical_solution1 << "\t" << "set1 = " << solution << "\t" << "error1 = " << error1 << "\t\t"
				<< "A2 = " << analytical_solution2 << "\t" << "set2 = " << "\t" << solution2 << "\t" << "error2 = " << error2 << endl;
		}
		
	}

	// ------------------------
	// Implicit Upwind Solution
	// ------------------------
	
	// Important Variables
	const double u = 1.75;
	const double courant = (u * dt) / dx;

	// Initialize Thomas algorithm vectors & vectors for solutions
	vector<double> f(N, 0.0);
	vector<double> fn(N, 0.0);
	vector<double> sub_diag(N, -courant);
	vector<double> main_diag(N, 1 + courant);
	vector<double> sup_diag(N, 0.0);
	cout << "***********************" << endl << "Implicit Solution Set 1" << endl << "***********************" << endl;
	
	// *********************
	// Set 1 implicit upwind
	// Set initial conditions
	for (int i = 0; i < N + 1; i++) {
		f[i] = set1_initial(i);
	}

	// Time step loop
	for (double t = 0; t < max_t + dt; t += dt) {
		vector<double> fin(N, 0.0);
		
		// Boundary conditions
		fin[0] = 0;
		fin[N-1] = 1;

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
			cout << "t = " << t << "\t" << "x = " << x_values[i] << "\t\t" << "solution = " << fn[i] << endl;
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
	cout << "***********************" << endl << "Implicit Solution Set 2" << endl << "***********************" << endl;
	// *********************
	// Set 2 implicit upwind
	// Set initial conditions
	for (int i = 0; i < N + 1; i++) {
		f[i] = set2_initial(i);
	}

	// Time step loop
	for (double t = 0; t < max_t + dt; t += dt) {
		vector<double> fin(N, 0.0);

		// Boundary conditions
		fin[0] = 0;
		fin[N - 1] = 0;

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
			cout << "t = " << t << "\t" << "x = " << x_values[i] << "\t\t" << "solution = " << fn[i] << endl;
		}

		f = fn;
	}

	// -----------------------------
	// Richtmyer Multi-Step Solution
	// -----------------------------
	
	// Vectors for solutions
	vector<double> r(N, 0);
	vector<double> r_half(N, 0);
	vector<double> r_new(N, 0);

	// Set 1 Richtmyer
	initial1(r, dx);
	cout << "***********************" << endl << "Richtmyer Solution Set 1" << endl << "***********************" << endl;
	// Time step loop
	for (double t = 0; t < max_t + dt; t += dt) {
		// Apply boundaries
		boundary1(r, N);

		// Predictor step
		for (int i = 0; i < N - 1; i++) {
			r_half[i] = 0.5 * (r[i] + r[i + 1]) - (courant/2) * (r[i + 1] - r[i]);
			cout << "t = " << t << "\t" << "x = " << x_values[i] << "\t" << "half-step = " << r_half[i] << endl;
		}

		r_half[0] = 0;
		r_half[N + 1] = 1;

		// Corrector step
		for (int i = 0; i < N - 1; i++) {
			r_new[i] = r[i] - courant * (r_half[i] - r_half[i - 1]);
			analytical_solution1 = analytical_1(x, t);
			//analytical_solution2 = analytical_2(x, i);
			error1 = r_new[i] - analytical_solution1;
			//error2 = r_new[i] - analytical_solution2;
			cout << "t = " << t << "\t" << "x = " << x_values[i] << "\t" << "corrector = " << r_new[i] << endl;
		}
		r = r_new;
		
	}
	
	// Set 2 Richtmyer
	initial2(r, dx);
	cout << "***********************" << endl << "Richtmyer Solution Set 2" << endl << "***********************" << endl;
	// Time step loop
	for (double t = 0; t < max_t + dt; t += dt) {
		// Apply boundary conditions
		boundary2(r, N);

		// Predictor step
		for (int i = 0; i < N - 1; i++) {
			r_half[i] = 0.5 * (r[i] + r[i + 1]) - (courant / 2) * (r[i + 1] - r[i]);
			cout << "t = " << t << "\t" << "x = " << x_values[i] << "\t" << "half-step = " << r_half[i] << endl;
		}

		r_half[0] = 0;
		r_half[N] = 0;

		//Corrector step
		for (int i = 0; i < N - 1; i++) {
			r_new[i] = r[i] - courant * (r_half[i] - r_half[i - 1]);
			analytical_solution2 = analytical_2(x, t);
			error2 = r_new[i] - analytical_solution2;
			cout << "t = " << t << "\t" << "x = " << x_values[i] << "\t" << "corrector = " << r_new[i] << endl;
		}
	}
	return 0;
}