/*
	testEn
*/

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_expint.h>
#include "simple.h"

using namespace std;

int main() {

cout << "testEn" << endl;

int n = -10, max = 1e8;
number x;

cout << "n = " << n << endl;
cout << "max = " << max << endl;

for (int j=0; j<max; j++) {
	x = 0.0001*(number)(j+1.0);
	gsl_sf_exp(x);
	//gsl_sf_expint_E1(x);
	//gsl_sf_expint_En(n,x);
}

return 0;
}
