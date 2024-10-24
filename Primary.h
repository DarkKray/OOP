#pragma once
#include <cmath>
#include <iostream>
#include <random>
#include <math.h>
#include <string>
#include <cstdlib>
#include "fstream"

using namespace std;

class Primary {
private:
	double mu = 0;
	double lambda = 1;
	double n = 1;
	double x = 0;

	unsigned long long factorial(const int &n);

	double calculateThirdMoment(const double& a, const double& b, const int& n_intervals,
		const double& mean);

	double calculateMean(const double& a, const double& b, const int& n_intervals);

public:
	Primary(int new_n, double mu, double lambda);
	void set_n(int &n);
	void set_mu(double &mu);
	void set_lambda(double &lambda);
	double get_n();
	double get_mu();
	double get_lambda();
	double Laplace_function(double x);
	double* get_metrics();
	double pseudo_number();
	double Laplace_random_value();
	void load(std::string &fileway);
	void save(string &fileway);
};
