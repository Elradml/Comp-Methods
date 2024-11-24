// Code to calculate the norms for all schemes for both set 1 and set 2
#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
using namespace std;

// Set 1 Norm calculations
void set1_norms(vector<double>& veupwind1, vector<double>& viupwind1, vector<double>& vlax1, vector<double>& vricht1, vector<double>& vanalytical1, double N) {
	// Norm 1 calculations
	// Necessary variables and vectors
	double norm1_e1 = 0;
	double norm1_i1 = 0;
	double norm1_l1 = 0;
	double norm1_r1 = 0;

	double norm2_e1 = 0;
	double norm2_i1 = 0;
	double norm2_l1 = 0;
	double norm2_r1 = 0;

	vector<double> eup_error1;
	vector<double> iup_error1;
	vector<double> lax_error1;
	vector<double> richt_error1;
	// Loop will calculate all norm 1, 2, and 3 values for each scheme
	for (int i = 0; i < N; i++) {
		// Explicit Upwind
		norm1_e1 += abs(veupwind1[i] - vanalytical1[i]);
		norm2_e1 += pow(veupwind1[i] - vanalytical1[i], 2);
		norm2_e1 = sqrt(norm2_e1);
		// Implicit Upwind
		norm1_i1 += abs(viupwind1[i] - vanalytical1[i]);
		norm2_i1 += pow(viupwind1[i] - vanalytical1[i], 2);
		norm2_i1 = sqrt(norm2_i1);
		
		// Lax Wendroff
		norm1_l1 += abs(vlax1[i] - vanalytical1[i]);
		norm2_l1 += pow(vlax1[i] - vanalytical1[i], 2);
		norm2_l1 = sqrt(norm2_l1);

		// Richtmyer
		norm1_r1 += abs(vricht1[i] - vanalytical1[i]);
		norm2_r1 += pow(vricht1[i] - vanalytical1[i], 2);
		norm2_r1 = sqrt(norm2_r1);

		// Error calculations for norm 3
		eup_error1.push_back(veupwind1[i] - vanalytical1[i]);
		iup_error1.push_back(viupwind1[i] - vanalytical1[i]);
		lax_error1.push_back(vlax1[i] - vanalytical1[i]);
		richt_error1.push_back(vricht1[i] - vanalytical1[i]);
	}
	// Norm 3 calculations
	auto norm3_e1 = max_element(eup_error1.begin(), eup_error1.end());
	auto norm3_i1 = max_element(iup_error1.begin(), iup_error1.end());
	auto norm3_l1 = max_element(lax_error1.begin(), lax_error1.end());
	auto norm3_r1 = max_element(richt_error1.begin(), richt_error1.end());

	// Norm outputs for set 1
	cout << "SET 1*****************" << endl << "Explicit Upwind:" << endl << "N1 = " << norm1_e1 << endl << "N2 = " << norm2_e1 << endl << "N3 = " << *norm3_e1 << endl; // Asterix next to norm 3 pointer to give value it points to
	cout << "Implciit Upwind:" << endl << "N1 = " << norm1_i1 << endl << "N2 = " << norm2_i1 << endl << "N3 = " << *norm3_i1 << endl;
	cout << "Lax-Wendroff:" << endl << "N1 = " << norm1_l1 << endl << "N2 = " << norm2_l1 << endl << "N3 = " << *norm3_l1 << endl;
	cout << "Richtmyer:" << endl << "N1 = " << norm1_r1 << endl << "N2 = " << norm2_r1 << endl << "N3 = " << *norm3_r1 << endl << endl;
}

// Set 2 Norm calculations
void set2_norms(vector<double>& veupwind2, vector<double>& viupwind2, vector<double>& vlax2, vector<double>& vricht2, vector<double>& vanalytical2, double N) {
	// Norm 2 calculations
	// Necessary variables
	double norm1_e2 = 0;
	double norm1_i2 = 0;
	double norm1_l2 = 0;
	double norm1_r2 = 0;

	double norm2_e2 = 0;
	double norm2_i2 = 0;
	double norm2_l2 = 0;
	double norm2_r2 = 0;

	vector<double> eup_error2;
	vector<double> iup_error2;
	vector<double> lax_error2;
	vector<double> richt_error2;
	// Loop will calculate all norm 1, 2, and 3 values for each scheme
	for (int i = 0; i < N; i++) {
		// Explicit Upwind
		norm1_e2 += abs(veupwind2[i] - vanalytical2[i]);
		norm2_e2 += pow(veupwind2[i] - vanalytical2[i], 2);
		norm2_e2 = sqrt(norm2_e2);

		// Implicit Upwind
		norm1_i2 += abs(viupwind2[i] - vanalytical2[i]);
		norm2_i2 += pow(viupwind2[i] - vanalytical2[i], 2);
		norm2_i2 = sqrt(norm2_i2);

		// Lax Wendroff
		norm1_l2 += abs(vlax2[i] - vanalytical2[i]);
		norm2_l2 += pow(vlax2[i] - vanalytical2[i], 2);
		norm2_l2 = sqrt(norm2_l2);

		// Richtmyer
		norm1_r2 += abs(vricht2[i] - vanalytical2[i]);
		norm2_r2 += pow(vricht2[i] - vanalytical2[i], 2);
		norm2_r2 = sqrt(norm2_r2);

		// Error calculations for norm 3
		eup_error2.push_back(veupwind2[i] - vanalytical2[i]);
		iup_error2.push_back(viupwind2[i] - vanalytical2[i]);
		lax_error2.push_back(vlax2[i] - vanalytical2[i]);
		richt_error2.push_back(vricht2[i] - vanalytical2[i]);
	}
	// Norm 3 calculations
	auto norm3_e2 = max_element(eup_error2.begin(), eup_error2.end());
	auto norm3_i2 = max_element(iup_error2.begin(), iup_error2.end());
	auto norm3_l2 = max_element(lax_error2.begin(), lax_error2.end());
	auto norm3_r2 = max_element(richt_error2.begin(), richt_error2.end());

	// Norm outputs for set 2
	cout << "SET 2***************" << endl << "Explicit Upwind:" << endl << "N1 = " << norm1_e2 << endl << "N2 = " << norm2_e2 << endl << "N3 = " << *norm3_e2 << endl;
	cout << "Implciit Upwind:" << endl << "N1 = " << norm1_i2 << endl << "N2 = " << norm2_i2 << endl << "N3 = " << *norm3_i2 << endl;
	cout << "Lax-Wendroff:" << endl << "N1 = " << norm1_l2 << endl << "N2 = " << norm2_l2 << endl << "N3 = " << *norm3_l2 << endl;
	cout << "Richtmyer:" << endl << "N1 = " << norm1_r2 << endl << "N2 = " << norm2_r2 << endl << "N3 = " << *norm3_r2 << endl;
}