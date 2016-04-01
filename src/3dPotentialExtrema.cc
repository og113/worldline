/*-------------------------------------------------------------------------------------------------------------------------
	3dPotentialExtrema
		program to investigate the extrema of the 3d potential, as in email from Arttu on 22/3/16
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <string>
#include "simple.h"
#include "gsl_extras.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - rootFinders
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/


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
