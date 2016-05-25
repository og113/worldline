/*
	highTemp
		calculates bounce worldline for charged scalar pair production at high temperature in nonrelativistic approximation
		does so by directly integrating equations of motion
*/

#include <ctime>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_integration.h>
#include "simple.h"
#include "parameters.h"
#include "genloop.h"
#include "nrloop.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - integrand
		1 - getting inputs from argv
		2 - defining basic quantities
		3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------------------------------------------------------
	0. integrand
----------------------------------------------------------------------------------------------------------------------------*/

struct paramsIntegrandStruct {
	double E;
	double kappa;
};

double Integrand (double x, void* parameters) {
	struct paramsIntegrandStruct* params = (struct paramsIntegrandStruct*)parameters;
	double E = params->E;
	double kappa = params->kappa;
	if ((E + 2.0 - x - kappa/4.0/PI/x)<0)
	cerr << "Integrand error: sqrt(<0)" << endl;
	return 1.0/sqrt(E + 2.0 - x - kappa/4.0/PI/x);
}

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

// data to print
string inputsFile = "inputs4";
bool verbose = true;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

cout << "using inputs file " << inputsFile << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	2 - defining basic quantities
-------------------------------------------------------------------------------------------------------------------------*/

#define dim 4

ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}

// parameter loops
uint Npl = pr.totalSteps();
if (verbose)
	cout << "looping over " << Npl << " steps" << endl;

// starting loop
for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pl>0)
		p = pr.position(pl);
		
	uint N = pow(2,p.K);

	if (verbose) {
		cout << "generating finite temperature bounce of ";
		cout << N << " points in " << dim << " dimensions ";
		cout << "using non-relativistic approximation" << endl;
	}

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
	string file;
	uint Seed = time(NULL);
	Loop<dim> loop(p.K,Seed);
	Metropolis<dim> met(loop,p,Seed);
	Point<dim> p0, dpz, dpt;
	number E = p.P4;
	number kappa = pow(p.G,3)*p.B;
	if ((kappa/4.0/PI)>pow(1.0-p.P4/2.0,2)) {
		cerr << "highTemp error: kappa(" << kappa << ") is too large" << endl;
		return 1;
	}
	number rL = (1.0-p.P4/2.0) - sqrt(pow((1.0-p.P4/2.0),2) - kappa/4.0/PI);
	number rR = (1.0-p.P4/2.0) + sqrt(pow((1.0-p.P4/2.0),2) - kappa/4.0/PI);
	if (verbose)
		cout << "rL = " << rL << ", rR = " << rR << endl;
	
	number beta, r, t;
	
	// setting up integration
	int workspace_size = 1000;
	number tolAbs = 0.0;
	number tolRel = 1.0e-7;
	number error;
	paramsIntegrandStruct params;
	params.E = E;
	params.kappa = kappa;
	gsl_function F;
	F.function = &Integrand;
	F.params = &params;
	number singularities[2];
	singularities[0] = rL;
	singularities[1] = rR;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (workspace_size);

	// finding beta
	gsl_integration_qagp(&F, singularities, 2, tolAbs, tolRel, workspace_size, w, &beta, &error);
	beta *= 2.0;
	if (verbose)
		cout << "E = " << E << ", beta = " << beta << ", error = " << error << endl;

	file = "data/highTemp/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(kappa)\
														+"_beta_"+nts(beta)+"_rank_"+nts(0)+".dat";
	
	for (uint k=0; k<N/4; k++) {
		r = rL + (rR-rL)*k/(number)(N/4.0-1.0);
		if (k>0) {
			singularities[1] = r;
			gsl_integration_qagp(&F, singularities, 2, tolAbs, tolRel, workspace_size, w, &t, &error);
		}
		
		dpz[dim-2] = r/2.0;
		dpt[dim-1] = -beta/2.0 + t;

		loop[k] = p0+dpz+dpt;
		loop[N/2-1-k] = p0+dpz-dpt;
		loop[k+N/2] = p0-dpz-dpt;
		loop[N-1-k] = p0-dpz+dpt;
	}
	
	cout << "printing to " << file << endl;
	loop.save(file);
	loop.clear();
	
	gsl_integration_workspace_free (w);

}

return 0;
}
