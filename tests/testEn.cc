/*
	testEn
*/

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#include "simple.h"

using namespace std;

int main() {

cout << "testEn" << endl;

int n = -100, max = 1e6;
number x, y = 0.0;

cout << "n = " << n << endl;
cout << "max = " << max << endl;

for (int j=0; j<max; j++) {
	x = 0.001*(number)(j+1.0);
	//y = gsl_sf_exp(x);
	//y = gsl_sf_expint_E1(x);
	//y = gsl_sf_expint_En(n,x);
	y = gsl_sf_gamma_inc(1.0-n,x)*pow(x,n-1.0);
	//cout << y << endl;
}
cout << y << endl;

return 0;
}
