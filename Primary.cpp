#include "Primary.h"
#include "string"
#include "fstream"

using namespace std;


Primary::Primary(int n, double mu, double lambda) {
	this->n = n;
	this->mu = mu;
	this->lambda = lambda;
}
void Primary::set_n(int &n_new) {
    this->n = n_new;

}
void Primary::set_mu(double &mu_new) {
    this->mu = mu_new;
}
void Primary::set_lambda(double &lambda_new) {
    this->lambda = lambda_new;
}
double Primary::get_n() {
	return this->n;
}
double Primary::get_mu() {
	return this->mu;
}
double Primary::get_lambda() {
	return this->lambda;
}

unsigned long long Primary::factorial(const int &n) {
    unsigned long long result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}
double Primary::Laplace_function(double x) {
	this->x = x;
    double value = exp(-abs(this->x - mu)) / (factorial(n - 1) * pow(2, n) * lambda);
    double sum = 0;
    for (int i = 0; i <= n - 1; i++) {
        sum += (factorial(n - 1 + i) / (factorial(n - 1 - i) * 
            factorial(i))) * (pow(abs(this->x - mu), n - 1 - i) / pow(2, i));
        
    }
    value = value * sum;
    return value;
}

double Primary::calculateThirdMoment(const double& a, const double& b, const int& n_intervals, 
    const double& mean) {

    double h = (b - a) / n_intervals;  // Шаг интегрирования
    double sum = 0.5 * (pow(a - mean, 3) * Laplace_function(a) + 
        pow(b - mean, 3) * Laplace_function(b));  // Начальные значения

    for (int i = 1; i < n_intervals; ++i) {
        double x = a + i * h;
        sum += pow(x - mean, 3) * Laplace_function(x);  // (x - mean)^3 * f(x)
    }

    return sum * h;
}

double Primary::calculateMean(const double& a, const double& b, const int& n_intervals) {
    double h = (b - a) / n_intervals;  // Шаг интегрирования
    double sum = 0.5 * (a * Laplace_function(a) + b * Laplace_function(b));  // Начальные значения

    for (int i = 1; i < n_intervals; ++i) {
        double x = a + i * h;
        sum += x * Laplace_function(x);  // x * f(x)
    }

    return sum * h;
}

double* Primary::get_metrics() {
    double n1 = n;
    double* result = new double[4];
    result[1] = 2 * n1; // дисперсия
    result[3] = 3 / n1; // эксцесс
    double a = -25;
    double b = 25;
    int n_intervals = 10000;

    double mean = calculateMean(a, b, n_intervals); // мат ожидание
    result[0] = mean;
    double thirdMoment = calculateThirdMoment(a, b, n_intervals, mean);
    double stddev = sqrt(result[1]);
    result[2] = thirdMoment / pow(stddev, 3); // коэф. ассиметрии

    return result;
}

double Primary::pseudo_number() {
    double random;
    do random = (double)rand() / RAND_MAX;
    while (random == 0. || random == 1.);
    return random;
}

double Primary::Laplace_random_value() {
    double p1 = 1;
    double p2 = 1;
    for (int i = 1; i <= n; i++) {
        double rand = pseudo_number();
        if (rand <= 0.5) {
            p1 = p1 * 2 * rand;
        }
        else {
            p2 = p2 * 2 * (1 - rand);
        }
    }
    return log(p1 / p2) * lambda + mu;
}

void Primary::load(string &fileway) {
    ifstream input(fileway);
    int tempNu;
    double tempMu, tempLambda;
    input >> tempNu;
    input >> tempMu;
    input >> tempLambda;
    input.close();  
    set_n(tempNu);
    set_mu(tempMu);
    set_lambda(tempLambda);
    cout << tempNu << n << endl;
    cout << tempMu << mu << endl;
    cout << tempLambda << lambda << endl;
    return;
}

void Primary::save(string& fileway) {
    ofstream output(fileway);
    output << n << "\n";
    output << mu << "\n";
    output << lambda << "\n";
    output.close();
    return;
}