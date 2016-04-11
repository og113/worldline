/*-------------------------------------------------------------------------------------------------------------------------
	declarations for 3d potential function and its derivatives
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __3DPOTENTIAL_H_INCLUDED__
#define __3DPOTENTIAL_H_INCLUDED__

#include <cmath>
#include "simple.h"
#include "parameters.h"

using namespace std;

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

struct params_for_V {
	double a;
	double kappa;
	void setFromParameters(const Parameters&);
};

/*-------------------------------------------------------------------------------------------------------------------------
	2. potential functions
		
	N.B. defined to work with gsl
-------------------------------------------------------------------------------------------------------------------------*/

double V0 (double y, void * parameters);

double dV0 (double y, void * parameters);

double ddV0 (double y, void * parameters);

double mV0 (double y, void * parameters);

double mdV0 (double y, void * parameters);

double mddV0 (double y, void * parameters);

double V4 (double y, void * parameters);

double dV4 (double y, void * parameters);

double ddV4 (double y, void * parameters);

double mV4 (double y, void * parameters);

double mdV4 (double y, void * parameters);

double mddV4 (double y, void * parameters);

#endif // __3DPOTENTIAL_H_INCLUDED__
