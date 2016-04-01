/*-------------------------------------------------------------------------------------------------------------------------
	definitions of functions giving potentials and their derivatives
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include "simple.h"
#include "parameters.h"
#include "3dpotential.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - parameters struct
	2 - potential functions
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. parameters struct
		
	N.B. defined to work with gsl
-------------------------------------------------------------------------------------------------------------------------*/

void params_for_V::setFromParameters(const Parameters& p) {
	a = p.Epsi;
	kappa = pow(p.G,3)*p.B;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. potential functions
		
	N.B. defined to work with gsl
-------------------------------------------------------------------------------------------------------------------------*/

double V (double y, void * parameters) {
	struct params_for_V* params = (struct params_for_V *)parameters;
	double y2 = y*y;
	double a = params->a;
	double a2 = pow(a,2);
	double kappa = params->kappa;
	return 2.0 - y - (kappa/4.0/PI)*( 1.0/sqrt(y2 + a2) - exp(-y2/a2)/a );
}

double dV (double y, void * parameters) {
	struct params_for_V* params = (struct params_for_V *)parameters;
	double y2 = y*y;
	double a = params->a;
	double a2 = pow(a,2);
	double kappa = params->kappa;
	return -1.0 - (kappa/4.0/PI)*( -y/pow(y2 + a2,1.5) + 2.0*y*exp(-y2/a2)/pow(a,3) );
}

double ddV (double y, void * parameters) {
	struct params_for_V* params = (struct params_for_V *)parameters;
	double y2 = y*y;
	double a = params->a;
	double a2 = pow(a,2);
	double kappa = params->kappa;
	return -(kappa/4.0/PI)*( (2.0*y2 - a2)/pow(y2 + a2,2.5) + (2.0/pow(a,3) - 4.0*y2/pow(a,5))*exp(-y2/a2) );
}

