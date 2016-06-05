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
#include <gsl/gsl_spline.h>
#include "simple.h"
#include "parameters.h"
#include "genloop.h"
#include "gsl_extras.h"
#include "nrloop.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - functions to use in main
		1 - getting inputs from argv
		2 - defining basic quantities
		3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------------------------------------------------------
	0. functions to use in main
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

struct paramsdS2Struct {
	paramsIntegralStruct paramsIntegral;
	double E;
	double L;
	uint M;
	double yi;
};

double Integrand (double x, void* parameters) {
	struct paramsIntegrandStruct* params = (struct paramsIntegrandStruct*)parameters;
	double E = params->E;
	double kappa = params->kappa;
	if ((E + 2.0 - x - kappa/4.0/PI/x)<0)
	cerr << "Integrand error: sqrt(<0)" << endl;
	return 2.0/sqrt(-E + 2.0 - x - kappa/4.0/PI/x);
}

double LIntegrand (double x, void* parameters) {
	struct paramsIntegrandStruct* params = (struct paramsIntegrandStruct*)parameters;
	double E = params->E;
	double kappa = params->kappa;
	if ((E + 2.0 - x - kappa/4.0/PI/x)<0)
	cerr << "LIntegrand error: sqrt(<0)" << endl;
	return sqrt(1.0 + 0.25*(2.0 - E - x - kappa/4.0/PI/x))/sqrt(2.0 - E - x - kappa/4.0/PI/x);
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

double dS2 (double y, void* parameters) {

	// getting parameters
	struct paramsdS2Struct* params = (struct paramsdS2Struct*)parameters;
	double E = params->E;
	double L = params->L;
	double yi = params->yi;
	uint M = params->M;
	struct paramsIntegralStruct paramsIntegral = params->paramsIntegral;
	paramsIntegral.a = yi;
	paramsIntegral.b = y;
	
	return pow(TIntegral(E,&paramsIntegral)/2.0,2) + pow((y-yi)/2.0,2) - pow(L/(number)M,2);
}

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

// data to print
string inputsFile = "inputs4";
bool verbose = true;
bool fixBeta = true;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixT")==0 || id.compare("fixBeta")==0 || id.compare("beta")==0) \
							fixBeta = (stn<uint>(argv[2*j+2])!=0);
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
	number beta = ((p.T)>sqrt(MIN_NUMBER) && fixBeta? 1.0/(p.T): 0.0);
	if (E>Ethreshold) {
		cerr << "highTemp error: E(" << E << ") above threshold(" << Ethreshold << ")" << endl;
		return 1;
	}
	number rL = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	number rR = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	if (verbose && fixBeta)
		cout << "beta = " << beta << endl;
	if (verbose && !fixBeta) 
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
		if (verbose)
			cout << "Left0 = " << BetaZeroIntegral(Emin,&params) << ", Right0 = " << BetaZeroIntegral(Emax,&params) << endl;
		E = brentRootFinder(&Beta_gsl,Eguess,Emin,Emax,1.0e-7);
		rL = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
		rR = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
		params.a = rL;
		params.b = rR;
		if (verbose) {
			cout << "rL = " << rL << ", rR = " << rR << endl;
			cout << "rHigh = " << sqrt(kappa/4.0/PI) << endl;
		}
	}
	else {
		beta = TIntegral(E,&params);
		params.beta = beta;
	}
	if (verbose)
		cout << "E = " << E << ", beta = " << beta << endl;

	file = "data/highTemp/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/highTemp_Kappa_"+nts(kappa)\
														+"_T_"+nts(1.0/beta)+"_rank_"+nts(0)+".dat";
	
	// calculating length
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(params.workspace_size);
	paramsIntegrandStruct paramsIntegrand;
	paramsIntegrand.E = E;
	paramsIntegrand.kappa = kappa;
	gsl_function FL;
	FL.function = &LIntegrand;
	FL.params = &paramsIntegrand;
	number L, errorL;
	gsl_integration_qags (&FL, rL, rR, params.tolAbs, params.tolRel, params.workspace_size, w, &L, &errorL); 
	gsl_integration_workspace_free (w);
	cout << "L = " << 4.0*L << endl;
	
	/*// getting interpolating function
	vector<number> r(N/4), t(N/4);	
	for (uint k=0; k<N/4; k++) {
		r[k] = rL + (rR-rL)*k/(number)(N/4.0-1.0);
		if (k>0) {
			params.b = r[k];
			t[k] = TIntegral(E,&params)/2.0;
		}
		else
			t[k] = 0.0;
	}	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N/4);
	gsl_spline_init (spline, &(t[0]), &(r[0]), N/4);*/
    
    // output
    number y = rL, ti = 0.0, t = ti;
    paramsdS2Struct paramsdS2;
    paramsdS2.paramsIntegral = params;
    paramsdS2.E = E;
    paramsdS2.L = L;
    paramsdS2.M = N/4;
    for (uint k=0; k<N/4; k++) {
    	/*ti = (beta/2.0)*(number)(k + 0.5)/(number)(N/4.0);
    	ri = gsl_spline_eval(spline, ti, acc); // use if want narrowest part at boundaries
    	//ri = gsl_spline_eval(spline, beta/2.0 - ti, acc); // use if want widest part at boundaries*/
    	ti = t;
    	if (k==0) {
    		paramsdS2.L = L/2.0;
    		paramsdS2.yi = rL;
    		gsl_function ds2_gsl;
			ds2_gsl.params = &paramsdS2;
			ds2_gsl.function = &dS2;
			number yMin = rL+params.tolRel;
			number yMax = rR;
			number yGuess = rL+(rR-rL)*2.0/(number)N;			
			y = brentRootFinder(&ds2_gsl,yGuess,yMin,yMax,params.tolRel);
    	}
    	else {
    		paramsdS2.L = L;
    		paramsdS2.yi = y;
    		gsl_function ds2_gsl;
			ds2_gsl.params = &paramsdS2;
			ds2_gsl.function = &dS2;
			number yMin = y+params.tolRel;
			number yMax = rR;
			number yGuess = y+(rR-rL)*4.0/(number)N;
			y = brentRootFinder(&ds2_gsl,yGuess,yMin,yMax,params.tolRel);
    	}
    	
    	params.a = rL;
    	params.b = y;
    	
    	t = TIntegral(E,&params);
    	
    	dpz[dim-2] = y/2.0;
		dpt[dim-1] = -beta/2.0 + t/2.0;

		loop[k] = p0+dpz+dpt;
		loop[N/2-1-k] = p0+dpz-dpt;
		loop[k+N/2] = p0-dpz-dpt;
		loop[N-1-k] = p0-dpz+dpt;
    }

    //gsl_spline_free(spline);
    //gsl_interp_accel_free (acc);

	cout << "printing to " << file << endl;
	loop.save(file);
	loop.clear();

}

return 0;
}
