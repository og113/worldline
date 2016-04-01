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

double V (double r, void * parameters) {
	struct params_for_V* params = (struct params_for_V *)parameters;
	double a = params->a;
	double kappa = params->kappa;
	return 2.0 - r - ((-(1.0/(a*exp(pow(r,2)/pow(a,2)))) + 1.0/sqrt(pow(a,2) + pow(r,2)))*kappa)/(4.0*PI);
}

double dV (double r, void * parameters) {
	struct params_for_V* params = (struct params_for_V *)parameters;
	double a = params->a;
	double kappa = params->kappa;
	return  - 1.0 - (((2.0*r)/(pow(a,3)*exp(pow(r,2)/pow(a,2))) - r/pow(pow(a,2) + pow(r,2),1.5))*kappa)/(4.0*PI);
}

double ddV (double r, void * parameters) {
	struct params_for_V* params = (struct params_for_V *)parameters;
	double a = params->a;
	double kappa = params->kappa;
	return  - ((2.0/(pow(a,3)*exp(pow(r,2)/pow(a,2))) \
			 - (4.0*pow(r,2))/(pow(a,5)*exp(pow(r,2)/pow(a,2))) \
			 + (-pow(a,2) + 2.0*pow(r,2))/pow(pow(a,2) \
			 + pow(r,2),2.5))*kappa)/(4.0*PI);
}

