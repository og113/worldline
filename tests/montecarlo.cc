// quick test to try to approximate a gaussian integral with monte carlo techniques

#include <ctime>
#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_randist.h> 	// Distributions of random numbers
#include "simple.h"

using namespace std;

int main() {

uint K = 2, N = 1000;
number sigma = 1.0;
number result = 0.0, x;
uint seed = time(NULL);
cout << "give power of moment: ";
cin >> K;
cout << "give number of terms: ";
cin >> N;

gsl_rng * generator;
generator = gsl_rng_alloc(gsl_rng_taus);
gsl_rng_set(generator,seed);

for (uint j=0; j<N; j++) {
	x = gsl_ran_gaussian_ziggurat(generator, sigma);
	result += pow(x,K);
}
result /= (number)N;

cout << "result = " << result << endl;

delete[] generator;

return 0;
}
