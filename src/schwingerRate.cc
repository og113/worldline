/*
	schwinger rate
		quick program to calculate schwinger rate for gnuplot to plot
*/

#include <cmath>
#include <iostream>
#include <string>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_dilog.h>
#include "simple.h"

using namespace std;

int main(int argc, char** argv) {

number x, argument, result;

if (argc!=2)  {
	cout << "Type argument, eE/m^2, to run schwingerRate: ";
	cin >> x;
}
else {
	x = stn<number>(argv[1]);
}

argument = -gsl_sf_exp(-pi/x);
result = pow(x,2)/(16.0*pow(pi,3));
result *= gsl_sf_dilog(argument);

cout << result << endl;

return 0;
}
