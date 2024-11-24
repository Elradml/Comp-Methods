/* A program which is intended to examine the application of numerical 
solutions to partial differencial equations using the linear advection problem
* using the explicit and implicit upwind, Lax-Wendroff, and Ritchmyer multi-step methods
using two sets of initial/boundary conditions */
#include <iostream>
#include <vector>
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
	r[N-1] = 1;
}

// Boundary conditions for Richtmyer Set 2
void boundary2(vector<double>& r, double N) {
	r[0] = 0;
	r[N-1] = 0;
}

void set1_norms(vector<double>& veupwind1, vector<double>& viupwind1, vector<double>& vlax1, vector<double>& vricht1, vector<double>& vanalytical1, double N);
void set2_norms(vector<double>& veupwind2, vector<double>& viupwind2, vector<double>& vlax2, vector<double>& vricht2, vector<double>& vanalytical2, double N);

//Main funciton
int main() {
	// ------------------------
	// Explicit Upwind Solution
	// ------------------------
	
	//defining vectors
	vector<double> x_values;
	vector<double> veupwind1;
	vector<double> veupwind2;
	vector<double> vanalytical1;
	vector<double> vanalytical2;
	vector<double> verror1;
	vector<double> verror2;
	vector<double> time_step;
	//defining domain
	const int min_x = -50;
	const int max_x = 50;
	// User enters desired grid spacing
	std::cout << "Please enter the number of space grid points desired: ";
	double N;
	cin >> N;
	std::cout << N << endl;
	//Set delta t and total time
	double dt = 0.1;
	const int max_t = 10;
	//other variables necessary for nested loop
	double solution;
	double solution2;
	double x;
	double dx;
	double analytical_solution1;
	double analytical_solution2;
	double error1;
	double error2;
	std::cout << "***********************" << endl << "Explicit Solutions" << endl << "***********************" << endl;

	//Time step loop
	for (double t = 0; t < max_t; t += dt) {
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

			error1 = solution - analytical_solution1;
			error2 = solution2 - analytical_solution2;

			if (fabs(t - 5) < 1e-6 || fabs(t - 10) < 1e-6) {
				veupwind1.push_back(solution);
				veupwind2.push_back(solution2);
				vanalytical1.push_back(analytical_solution1);
				vanalytical2.push_back(analytical_solution2);
				verror1.push_back(error1);
				verror2.push_back(error2);
				time_step.push_back(t);
			}
		}
	}
	// Output values to .csv
	ofstream file("eupwind.csv");
	// Headers
	file << "Time step, X values, Set 1 Solution, Set 1 Analytical, Set 1 Error, Set 2 Solution, Set 2 Analytical, Set 2 Error\n";
	// Write data to file
	for (size_t i = 0; i < time_step.size(); i++) {
		file << time_step[i] << ", " << x_values[i] << "," << veupwind1[i] << "," << vanalytical1[i] << "," << verror1[i] << "," << veupwind2[i] << "," << vanalytical2[i] << "," << verror2[i] << "\n";
	}
	file.close();

	// ------------------------
	// Implicit Upwind Solution
	// ------------------------
	
	// Important Variables
	const double u = 1.75;
	const double courant = (u * dt) / dx;
	double temp;

	// Initialize Thomas algorithm vectors & vectors for solutions
	vector<double> f(N, 0.0);
	vector<double> fn(N, 0.0);
	vector<double> sub_diag(N, -courant);
	vector<double> main_diag(N, 1 + courant);
	vector<double> sup_diag(N, 0.0);
	vector<double> viupwind1;
	vector<double> viupwind2;
	vector<double> vierror1;
	std::cout << "***********************" << endl << "Implicit Solution Set 1" << endl << "***********************" << endl;
	// *********************
	// Set 1 implicit upwind
	// Set initial conditions
	
	for (int i = 0; i < N; i++) {
		temp = x_values[i];
		f[i] = set1_initial(temp);
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
		}

		// Correcting vector values for output
		if (fabs(t - 5) < 1e-6 || fabs(t - 10) < 1e-6) {
			for (int i = 0; i < N; i++) {
				error1 = fn[i] - vanalytical1[i];
				viupwind1.push_back(fn[i]);
				vierror1.push_back(error1);
			}
		}
		f = fn;
		

	}
	std::cout << "***********************" << endl << "Implicit Solution Set 2" << endl << "***********************" << endl;
	// Clearing vectors for set 2
	/*f.clear();
	fn.clear();
	sub_diag.clear();
	main_diag.clear();
	sup_diag.clear();*/
	
	// *********************
	// Set 2 implicit upwind
	// Set initial conditions
	for (int i = 0; i < N; i++) {
		temp = x_values[i];
		f[i] = set2_initial(temp);
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
		}
		if (fabs(t - 5) < 1e-6 || fabs(t - 10) < 1e-6) {
			for (int i = 0; i < N; i++) {
				error2 = fn[i] - vanalytical2[i];
				viupwind2.push_back(fn[i]);
				verror2.at(i) = error2;
			}
		}

		f = fn;
	}
	// Output values to .csv
	ofstream file1("iupwind.csv");
	// Headers
	file1 << "Time step, X values, Set 1 Solution, Set 1 Analytical, Set 1 Error, Set 2 Solution, Set 2 Analytical, Set 2 Error\n";
	// Write data to file
	for (size_t i = 0; i < time_step.size(); i++) {
		file1 << time_step[i] << ", " << x_values[i] << "," << viupwind1[i] << "," << vanalytical1[i] << "," << vierror1[i] << "," << viupwind2[i] << "," << vanalytical2[i] << "," << verror2[i] << "\n";
	}
	file1.close();

	std::cout << "Lax wendroff start" << endl;
	// ------------------------
	// Lax-Wendroff Solution
	// ------------------------
	vector<double> LaxWendroff_set1(N, 0.0); //Solution vector for Lax-Wendroff
	vector<double> LaxWendroff_set1_new(N, 0.0); //Temporary vector for updates

	vector<double> vlax1;
	vector<double> vlax2;
	vector<double> lerror1;
	vector<double> lerror2;

	//Initialise the soluton with intial conditions
	for (int i = 0; i < N; i++) {
		x = min_x + (max_x - min_x) * (i / (N - 1));
		LaxWendroff_set1[i] = set1_initial(x);
	}

	//Time step loop
	for (double t = 0; t < max_t; t += dt) {
		// Apply the Lax-Wendroff for the interior points
		for (int i = 1; i < N - 1; i++) {
			LaxWendroff_set1_new[i] = LaxWendroff_set1[i] - 0.5 * courant * (LaxWendroff_set1[i + 1] - LaxWendroff_set1[i - 1])
				+ 0.5 * courant * courant * (LaxWendroff_set1[i + 1] - 2 * LaxWendroff_set1[i] + LaxWendroff_set1[i - 1]);

			if (fabs(t - 5) < 1e-6 || fabs(t - 10) < 1e-6) {
				for (int i = 0; i < N; i++) {
					error1 = LaxWendroff_set1_new[i] - vanalytical1[i];
					vlax1.push_back(LaxWendroff_set1_new[i]);
					lerror1.push_back(error1);
				}
			}
		}

		//Boundary Conditions 
		LaxWendroff_set1_new[0] = 0;     //Left boundary
		LaxWendroff_set1_new[N - 1] = 1; //Right boundary

		//Upadte the solution
		LaxWendroff_set1 = LaxWendroff_set1_new;
		// Calculation of errors and setting up vectors for output set 1

	}


	// ------------------------
	// Lax-Wendroff Solution for Set2
	// ------------------------
	vector<double> LaxWendroff_set2(N, 0.0); // Solution vector for Lax-Wendroff
	vector<double> LaxWendroff_set2_new(N, 0.0); // Temporary vector for updates

	// Initialise the solution with initial conditions
	for (int i = 0; i < N; i++) {
		x = min_x + (max_x - min_x) * (i / (N - 1));
		LaxWendroff_set2[i] = set2_initial(x);
	}

	// Time step loop
	for (double t = 0; t < max_t; t += dt) {
		// Apply the Lax-Wendroff for the interior points
		for (int i = 1; i < N - 1; i++) {
			LaxWendroff_set2_new[i] = LaxWendroff_set2[i]
				- 0.5 * courant * (LaxWendroff_set2[i + 1] - LaxWendroff_set2[i - 1])
				+ 0.5 * courant * courant * (LaxWendroff_set2[i + 1] - 2 * LaxWendroff_set2[i] + LaxWendroff_set2[i - 1]);

			// Calculation of errors and setting up vectors for output set 2
			if (fabs(t - 5) < 1e-6 || fabs(t - 10) < 1e-6) {
				for (int i = 0; i < N; i++) {
					error2 = LaxWendroff_set2_new[i] - vanalytical2[i];
					vlax2.push_back(LaxWendroff_set2_new[i]);
					lerror2.push_back(error2);
				}
			}
		}

		// Boundary Conditions 
		LaxWendroff_set2_new[0] = 0;     // Left boundary
		LaxWendroff_set2_new[N - 1] = 1; // Right boundary

		// Update the solution
		LaxWendroff_set2 = LaxWendroff_set2_new;
	}


	// Output values to .csv
	ofstream file2("lax.csv");
	// Headers
	file2 << "Time step, X values, Set 1 Solution, Set 1 Analytical, Set 1 Error, Set 2 Solution, Set 2 Analytical, Set 2 Error\n";
	// Write data to file
	for (size_t i = 0; i < time_step.size(); i++) {
		file2 << time_step[i] << ", " << x_values[i] << "," << vlax1[i] << "," << vanalytical1[i] << "," << lerror1[i] << "," << vlax2[i] << "," << vanalytical2[i] << "," << lerror2[i] << "\n";
	}
	file2.close();

	// -----------------------------
	// Richtmyer Multi-Step Solution
	// -----------------------------
	
	// Vectors for solutions
	vector<double> r(N, 0);
	vector<double> r_half(N, 0);
	vector<double> r_new(N, 0);

	vector<double> vricht1;
	vector<double> vricht2;
	vector<double> vrerror1;
	vector<double> vrerror2;

	// Set 1 Richtmyer
	initial1(r, dx);
	std::cout << "***********************" << endl << "Richtmyer Solution Set 1" << endl << "***********************" << endl;
	// Time step loop
	for (double t = 0; t < max_t + dt; t += dt) {
		// Apply boundaries
		boundary1(r, N);
		// Predictor step
		for (int i = 0; i < N - 1; i++) {
			r_half[i] = 0.5 * (r[i] + r[i + 1]) - (courant/2) * (r[i + 1] - r[i]);
		}

		r_half[0] = 0;
		r_half[N - 1] = 1;
		// Corrector step
		for (int i = 1; i < N; i++) {
			r_new[i] = r[i] - courant * (r_half[i] - r_half[i - 1]);
		}
		// Vectors for output
		if (fabs(t - 5) < 1e-6 || fabs(t - 10) < 1e-6) {
			for (int i = 0; i < N; i++) {
				vricht1.push_back(r_new[i]);
				error1 = vricht1[i] - vanalytical1[i];
				vrerror1.push_back(error1);
			}
		}

		r = r_new;
		
	}
	
	// Set 2 Richtmyer
	initial2(r, dx);
	std::cout << "***********************" << endl << "Richtmyer Solution Set 2" << endl << "***********************" << endl;
	// Time step loop
	for (double t = 0; t < max_t + dt; t += dt) {
		// Apply boundary conditions
		boundary2(r, N);

		// Predictor step
		for (int i = 0; i < N - 1; i++) {
			r_half[i] = 0.5 * (r[i] + r[i + 1]) - (courant / 2) * (r[i + 1] - r[i]);
			//cout << "t = " << t << "\t" << "x = " << x_values[i] << "\t" << "half-step = " << r_half[i] << endl;
		}

		r_half[0] = 0;
		r_half[N - 1] = 0;

		//Corrector step
		for (int i = 1; i < N; i++) {
			r_new[i] = r[i] - courant * (r_half[i] - r_half[i - 1]);
		}

		// Vectors for output
		if (fabs(t - 5) < 1e-6 || fabs(t - 10) < 1e-6) {
			for (int i = 0; i < N; i++) {
				vricht2.push_back(r_new[i]);
				error2 = vricht2[i] - vanalytical2[i];
				vrerror2.push_back(error2);
			}
		}
	}
	// Output values to .csv
	ofstream file3("richt.csv");
	// Headers
	file3 << "Time step, X values, Set 1 Solution, Set 1 Analytical, Set 1 Error, Set 2 Solution, Set 2 Analytical, Set 2 Error\n";

	// Write data to file
	for (size_t i = 0; i < time_step.size(); i++) {

		file3 << time_step[i] << ", " << x_values[i] << "," << vricht1[i] << "," << vanalytical1[i] << "," << vrerror1[i] << "," << vricht2[i] << "," << vanalytical2[i] << "," << vrerror2[i] << "\n";
	}
	file3.close();

	// Call functions to calculate norms for both set 1 and 2
	set1_norms(veupwind1, viupwind1, vlax1, vricht1, vanalytical1, N);
	set2_norms(veupwind2, viupwind2, vlax2, vricht2, vanalytical2, N);
	return 0;
}