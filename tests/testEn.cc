/*
	testEn
*/

//#include <gsl/gsl_math.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#include "simple.h"

using namespace std;

int main() {

cout << "testEn" << endl;

int n = -100;
long int max = 1e8;
number x = 1.0, y = 0.0;

cout << "n = " << n << endl;
cout << "max = " << max << endl;

for (long int j=0; j<max; j++) {
	x = 1.0e-8*(number)(j+1.0);
	//y = exp(-x);
	y = cos(x);
	//y = gsl_sf_expint_E1(x);
	//y = gsl_sf_expint_En(n,x);
	//y = gsl_sf_gamma_inc(1.0-n,x)*pow(x,n-1.0);
	//cout << y << endl;
}
cout << y << endl;

return 0;
}
