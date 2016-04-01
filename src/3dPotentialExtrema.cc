/*-------------------------------------------------------------------------------------------------------------------------
	3dPotentialExtrema
		program to investigate the extrema of the 3d potential, as in email from Arttu on 22/3/16
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <string>
#include "simple.h"
#include "gsl_extras.h"
#include "3dpotential.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - loading inputs and argv
	2 - finding minima and maxima
	3 - printing results
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/


using namespace std;

int main(int argc, char** argv) {
/*-------------------------------------------------------------------------------------------------------------------------
	1 - loading inputs and argv
-------------------------------------------------------------------------------------------------------------------------*/

// argv options
bool verbose = true;
string inputsFile = "inputs4";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

if (verbose)
	cout << "using inputs file " << inputsFile << endl;

// loading inputs
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;
params_for_V params;
params.setFromParameters(p);

/*-------------------------------------------------------------------------------------------------------------------------
	2 - finding minima and maxima
-------------------------------------------------------------------------------------------------------------------------*/

// initial guesses and ranges
double max_guess = sqrt(params.kappa/4.0/PI);
double max_l = max_guess - 2.0*params.a;
double max_r = max_guess + 2.0*params.a;
double min_guess = params.a;
double min_l = 0.0;
double min_r = max_l;

// asigning functions
gsl_function dV_gsl;
dV_gsl.function = &dV;
dV_gsl.params = &params;

// finding extrema
double max = brentRootFinder(&dV_gsl,max_guess,max_l,max_r);
double min = brentRootFinder(&dV_gsl,min_guess,min_l,min_r);

// finding value of V at min and max
double Vmax = V(max,&params);
double Vmin = V(min,&params);

/*-------------------------------------------------------------------------------------------------------------------------
	3 - printing results	
-------------------------------------------------------------------------------------------------------------------------*/

printf("\n");
printf("%20s%20s%20s%20s%20s%20s\n","kappa","a","max","min","Vmax","Vmin");
printf("%20.16g%20.16g%20.16g%20.16g%20.16g%20.16g\n",params.kappa,params.a,max,min,Vmax,Vmin);
printf("\n");

return 0;
}
