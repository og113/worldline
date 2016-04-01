/*
	test3dpotential
*/

#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "gsl_extras.h"
#include "3dpotential.h"
#include "parameters.h"

using namespace std;

int main() {

cout << "test3dpotential:" << endl;

// inputs file
string pfile = "inputs4";
cout << "loading parameters from " << pfile << endl;
ParametersRange pr;
pr.load(pfile);
Parameters p = pr.Min;
params_for_V params;
params.setFromParameters(p);

cout << "a = " << params.a << endl;
cout << "kappa = " << params.kappa << endl;

// quick loop over y
double y = 0.0;
cout << left << setprecision(16);
cout << setw(25) << "y" << setw(25) << "V(y)" << setw(25) << "dV(y)" << setw(25) << "ddV(y)" << endl;
for (uint j=0; j<10; j++) {
	cout << setw(25) << y << setw(25) << V(y,&params) << setw(25) << dV(y,&params) << setw(25) << ddV(y,&params) << endl;
	y=0.01*(double)(j+1.0);
}

return 0;
}
