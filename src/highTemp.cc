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
#include "gsl_extras.h"
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

struct paramsIntegralStruct {
	double beta;
	double kappa;
	double a;
	double b;
	int workspace_size;
	number tolAbs;
	number tolRel;
};

double Integrand (double x, void* parameters) {
	struct paramsIntegrandStruct* params = (struct paramsIntegrandStruct*)parameters;
	double E = params->E;
	double kappa = params->kappa;
	if ((E + 2.0 - x - kappa/4.0/PI/x)<0)
	cerr << "Integrand error: sqrt(<0)" << endl;
	return 2.0/sqrt(-E + 2.0 - x - kappa/4.0/PI/x);
}

double TIntegral (double E, void* parameters) {
	// getting parameters
	struct paramsIntegralStruct* params = (struct paramsIntegralStruct*)parameters;
	double kappa = params->kappa;
	double a = params->a;
	double b = params->b;
	int workspace_size = params->workspace_size;
	number tolAbs  = params->tolAbs;
	number tolRel  = params->tolRel;
	
	// calculating other parameters
	paramsIntegrandStruct paramsIntegrand;
	paramsIntegrand.E = E;
	paramsIntegrand.kappa = kappa;
	gsl_function F;
	F.function = &Integrand;
	F.params = &paramsIntegrand;
	number singularities[2];
	singularities[0] = a;
	singularities[1] = b;
	
	// initializing workspace
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(workspace_size);
	
	// result
	number beta, error;

	// finding beta
	gsl_integration_qagp(&F, singularities, 2, tolAbs, tolRel, workspace_size, w, &beta, &error);
	
	// clearing workspace
	gsl_integration_workspace_free (w);
	
	return beta;
}

double BetaZeroIntegral (double E, void* parameters) {

	// getting parameters
	struct paramsIntegralStruct* params = (struct paramsIntegralStruct*)parameters;
	
	// fixing endpoints
	(*params).a = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - (params->kappa)/4.0/PI);
	(*params).b = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - (params->kappa)/4.0/PI);
	
	return (TIntegral(E,params) - params->beta);
}

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

// data to print
string inputsFile = "inputs4";
bool verbose = true;
bool fixBeta = false;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixBeta")==0 || id.compare("beta")==0) fixBeta = (stn<uint>(argv[2*j+2])!=0);
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
	number kappa = pow(p.G,3)*p.B;
	number Ethreshold = 2.0*(1.0-sqrt(kappa/4.0/PI));
	number E = (fixBeta? 1.0: p.P4);
	number beta = ((p.T*p.G*p.B)>sqrt(MIN_NUMBER) && fixBeta? 1.0/(p.T*p.G*p.B): 0.0);
	number r, t;
	if (E>Ethreshold) {
		cerr << "highTemp error: E(" << E << ") above threshold(" << Ethreshold << ")" << endl;
		return 1;
	}
	number rL = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	number rR = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	if (verbose) 
		cout << "rL = " << rL << ", rR = " << rR << endl;
	
	
	
	// params for integration
	paramsIntegralStruct params;
	params.beta = beta;
	params.kappa = kappa;
	params.a = rL;
	params.b = rR;
	params.workspace_size = 1e5;
	params.tolAbs = 0.0;
	params.tolRel = 1.0e-7;

	// finding beta or E
	if (fixBeta) {
		gsl_function Beta_gsl;
		Beta_gsl.params = &params;
		Beta_gsl.function = &BetaZeroIntegral;
		number Emin = 1.0e-7;
		number Emax = Ethreshold-1.0e-7;
		number Eguess = (Emax+Emin)/2.0;
		E = brentRootFinder(&Beta_gsl,Eguess,Emin,Emax,1.0e-7);
		rL = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
		rR = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	}
	else {
		beta = TIntegral(E,&params);
	}
	if (verbose)
		cout << "E = " << E << ", beta = " << beta << endl;

	file = "data/highTemp/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/highTemp_kappa_"+nts(kappa)\
														+"_T_"+nts(1.0/beta)+"_rank_"+nts(0)+".dat";
	
	for (uint k=0; k<N/4; k++) {
		r = rL + (rR-rL)*k/(number)(N/4.0-1.0);
		if (k>0) {
			params.b = rR;
			beta = TIntegral(E,&params);
		}
		else
			t = 0.0;
		
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

}

return 0;
}
