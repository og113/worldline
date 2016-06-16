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

struct paramsTotalStruct {
	double beta;
	double kappa;
	double a;
	double b;
	double E;
	double L;
	double f;
	int workspace_size;
	number tolAbs;
	number tolRel;
};

double TIntegrand (double x, void* parameters) {
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	double E = params.E;
	double kappa = params.kappa;
	if ((-E + 2.0 - x - kappa/4.0/PI/x)<0) {
		cerr << "TIntegrand error: sqrt(" << (-E + 2.0 - x - kappa/4.0/PI/x) << ")" << endl;
		cerr << "x = "<<  x << endl;
	}
	return 2.0/sqrt(-E + 2.0 - x - kappa/4.0/PI/x);
}

double TIntegral (double E, void* parameters) {
	// getting parameters
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	params.E = E;
	number a = params.a;
	number b = params.b;
	int workspace_size = params.workspace_size;
	number tolAbs  = params.tolAbs;
	number tolRel  = params.tolRel;
	
	if (abs(b-a)<MIN_NUMBER)
		return 0.0;
	
	// calculating other parameters
	gsl_function F;
	F.function = &TIntegrand;
	F.params = &params;
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
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	number kappa = params.kappa;
	params.E = E;
	
	// fixing endpoints
	params.a = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	params.b = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	
	return (TIntegral(E,&params) - params.beta);
}

double LIntegrand (double x, void* parameters) {
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	double E = params.E;
	double kappa = params.kappa;
	if ((-E + 2.0 - x - kappa/4.0/PI/x)<0){
		cerr << "LIntegrand error: sqrt(" << (-E + 2.0 - x - kappa/4.0/PI/x) << ")" << endl;
		cerr << "x = " << x << endl;
	}
	return sqrt(1.0 + 0.25*(2.0 - E - x - kappa/4.0/PI/x))/sqrt(2.0 - E - x - kappa/4.0/PI/x);
}

double LIntegral (double x, void* parameters) {
	// getting parameters
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	double a = params.a;
	double b = x;
	params.b = b;
	int workspace_size = params.workspace_size;
	number tolAbs  = params.tolAbs;
	number tolRel  = params.tolRel;
	
	if (abs(b-a)<MIN_NUMBER)
		return 0.0;
	
	// getting other parameters
	gsl_function F;
	F.function = &LIntegrand;
	F.params = &params;
	number singularities[2];
	singularities[0] = a;
	singularities[1] = b;
	
	// initializing workspace
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(workspace_size);
	
	// result
	number L, error;

	// finding beta
	gsl_integration_qagp(&F, singularities, 2, tolAbs, tolRel, workspace_size, w, &L, &error);
	
	// clearing workspace
	gsl_integration_workspace_free (w);
	
	return L;
}

double dS2 (double y, void* parameters) {

	// getting parameters
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	double E = params.E;
	double L = params.L;
	double yi = params.a;
	double f = params.f;
	params.b = y;
	
	return pow(TIntegral(E,&params)/2.0,2) + pow((y-yi)/2.0,2) - pow(L*f,2);
}

double dt2 (double y, void* parameters) {

	// getting parameters
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	double E = params.E;
	double L = params.L;
	double f = params.f;
	params.b = y;
	
	return pow(TIntegral(E,&params)/2.0,2) - pow(L*f,2);
}

double Lfway (double y, void* parameters) {

	// getting parameters
	struct paramsTotalStruct params = *((struct paramsTotalStruct*)parameters);
	double L = params.L;
	double f = params.f;
	
	return LIntegral(y,&params) - f*L;
}

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

// data to print
string inputsFile = "inputs4";
bool verbose = true;
bool fixBeta = true;
bool fixDS = true;
bool printInputs = false;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixT")==0 || id.compare("fixBeta")==0 || id.compare("beta")==0) \
							fixBeta = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixDS")==0 || id.compare("fixDs")==0 || id.compare("fixds")==0) \
							fixDS = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("printInputs")==0) printInputs = (stn<uint>(argv[2*j+2])!=0);
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

ParametersRange pr, prOut;
pr.load(inputsFile);
prOut = pr;
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
		if (fixDS)
			cout << "fixing ds^2" << endl;
	}

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
	
	// filename
	string file;
	
	// loop initialisation
	uint Seed = time(NULL);
	Loop<dim> loop(p.K,Seed);
	Metropolis<dim> met(loop,p,Seed);
	
	// basic parameters
	number kappa = pow(p.G,3)*p.B;
	if (verbose)
		cout << "kappa = " << kappa << endl;
	number Ethreshold = 2.0*(1.0-sqrt(kappa/4.0/PI));
	number E = (fixBeta? 1.0: p.P4);
	number beta = ((p.T)>sqrt(MIN_NUMBER) && fixBeta? 1.0/(p.T): 0.0);
	number betaMin = sqrt(2.0)*pow(PI,3.0/4.0)*pow(kappa,0.25);
	if (E>Ethreshold) {
		cerr << "highTemp error: E(" << E << ") above threshold(" << Ethreshold << ")" << endl;
		return 1;
	}
	if (verbose && fixBeta)
		cout << "beta = " << beta << endl;
	
	// endpoints (unchanged if !fixDS)
	number rL = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	number rR = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
	if (verbose && !fixBeta) 
		cout << "rL = " << rL << ", rR = " << rR << endl;
		
	// tolerances
	number tol1 = 1.0e-6, tol2 = 1.0e-7, tol3 = 1.0e-10, tol4 = 1.0e-12;
	
	// params for integration
	paramsTotalStruct params;
	params.beta = beta;
	params.kappa = kappa;
	params.a = rL;
	params.b = rR;
	params.workspace_size = 1e5;
	params.tolAbs = 0.0;
	params.tolRel = tol2;

	// finding beta or E
	if (fixBeta) {
		gsl_function Beta_gsl;
		Beta_gsl.params = &params;
		Beta_gsl.function = &BetaZeroIntegral;
		number Emin = tol1;
		number Emax = Ethreshold - tol1;
		number Eguess = (Emax+Emin)/2.0;
		E = brentRootFinder(&Beta_gsl,Eguess,Emin,Emax,tol4);
		rL = (1.0-E/2.0) - sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
		rR = (1.0-E/2.0) + sqrt(pow((1.0-E/2.0),2) - kappa/4.0/PI);
		params.E = E;
		params.a = rL;
		params.b = rR;
	}
	else {
		beta = TIntegral(E,&params);
		params.beta = beta;
	}
	if (verbose) {
		cout << "E = " << E << ", beta = " << beta << endl;
		cout << "Ethreshold = " << Ethreshold << ", betaMin = " << betaMin << endl;
		cout << "Ethreshold-E = " << Ethreshold-E << ", (pow(kappa,0.25)/1.11171)*(beta-betaMin) = " << (pow(kappa,0.25)/2.22341)*(beta-betaMin) << endl;
	}

	file = "data/highTemp/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/highTemp_Kappa_"+nts(kappa)\
														+"_T_"+nts(1.0/beta)+"_rank_"+nts(0)+".dat";
	
	// calculating length
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(params.workspace_size);
	gsl_function FL;
	FL.function = &LIntegrand;
	FL.params = &params;
	number L, errorL;
	gsl_integration_qags (&FL, rL, rR, params.tolAbs, tol3, params.workspace_size, w, &L, &errorL); 
	gsl_integration_workspace_free (w);
	params.L = L;
	cout << "L = " << 4.0*L << endl;
	if (verbose) {
		cout << "rL = " << rL << ", rR = " << rR << endl;
		cout << "need: dx << a << rL" << endl << "      (rR-rL)/2.0 << beta/2.0 < 1" << endl;
		cout << "      (rR-rL) << (rR+rL)/2.0" << endl;
		cout << "      (Ethreshold-E) << beta^2/((rR+rL)/2.0)" << endl;
		cout << "      (Ethreshold-E) << (rR+rL)/2.0)/4.0" << endl;
		cout << "have: " << 4.0*L/(number)N << " << " << p.Epsi << " << " << rL;
		cout << endl << "      " << (rR-rL)/2.0 << " << " << beta/2.0 << " < 1" << endl;
		cout << "      " << (rR-rL) << " << " << (rR+rL)/2.0 << endl;
		cout << "      " << (Ethreshold-E) << " << " << pow(beta,2)/((rR+rL)/2.0) << endl;
		cout << "      " << (Ethreshold-E) << " << " << (rR+rL)/2.0/4.0 << endl;
	}
	
	// printing suitable parameter range for given beta
	if (fixBeta && printInputs) {
		string inputsFileOut = "temp/inputsHighTemp_run_" + nts(pl);
		uint steps = 10;
		number muchMore = 3.0;
		number muchmuchMore = 10.0;
		number kappaMin = 1.0/(4.0*pow(PI,3)*pow(beta,4));
		prOut.Min = p;
		prOut.Max = p;
		fill((prOut.Steps).begin(),(prOut.Steps).end(),0.0);
		(prOut.Steps)[7] = steps;
		(prOut.Steps)[9] = steps;
		(prOut.Min).B = 1.2*kappaMin;
		(prOut.Max).B = muchmuchMore*1.2*kappaMin; // kappa max is difficult to calculate so just fudging it
		(prOut.Min).Epsi = muchMore*2.0*beta/(number)N;
		(prOut.Max).Epsi = sqrt((prOut.Max).B/4.0/PI)/muchMore;
		prOut.save(inputsFileOut);
		cout << "saved inputs to: " << inputsFileOut << endl;
	}	
    
    // output  
    Point<dim> p0, dpz, dpt;
    if (!fixDS) { //////// if not fixing ds2
    	
    	// getting interpolating function
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
		gsl_spline_init (spline, &(t[0]), &(r[0]), N/4);
		cout << "generated interpolating function" << endl;

    	number zi = rL, ti = 0.0;
    	for (uint k=0; k<N/4; k++) {
    	
			ti = (beta/2.0)*((number)k + 0.5)/((number)N/4.0);
			zi = 0.5*gsl_spline_eval(spline, ti, acc); // use if want narrowest part at boundaries
			//y = gsl_spline_eval(spline, beta/2.0 - ti, acc); // use if want widest part at boundaries
		
			params.a = rL;
			params.b = 2.0*zi;
			
			
			dpz[dim-2] = zi;
			dpt[dim-1] = -beta/2.0 + (beta/2.0)*((number)k + 0.5)/((number)N/4.0);

			loop[k] = p0+dpz+dpt;
			loop[N/2-1-k] = p0+dpz-dpt;
			loop[k+N/2] = p0-dpz-dpt;
			loop[N-1-k] = p0-dpz+dpt;
		}
		       
		gsl_spline_free(spline);
		gsl_interp_accel_free (acc);
    }
    else { //////// output if fixing ds2
    
    	// sorting params
    	params.L = L;
    	params.f = 0.5;
    
    	// integers
    	uint M = N/2;
    	uint J = p.K-1;
    	uint u, v, w;
    	
    	// numbers
    	number zi = 0.0, ti = 0.0;
    	
    	// for root finding
    	gsl_function Lfway_gsl;
		Lfway_gsl.params = &params;
		Lfway_gsl.function = &Lfway;
		number rMin, rMax, rGuess;
    
    	// first step
    	dpt[dim-1] = -beta/2.0;
    	dpz[dim-2] = rL/2.0;  // so narrowest part at boundaries
    	loop[0] = p0 + dpt + dpz;
    	
    	// second step
    	dpz[dim-2] = rR/2.0;
    	loop[M/2] = p0 + dpz;
    	
    	// following steps
    	for (uint r=2; r<=J; r++) {
    		for (uint s=0; s<pow(2,r-1); s++) {
    			u = (M/pow(2,r))*2*s;
    			v = (M/pow(2,r))*(2*s+1);
    			w = (M/pow(2,r))*(2*s+2);
    			
    			if (r==J)
    				u += 0;
    			
    			// loop[v] = halfway between loop[u] and pM or between loop[u] and loop[w]
    			if (v>M/2) {
    				params.a = (s==(pow(2,r-1)-1)? rL: 2.0*(loop[w])[dim-2]);
    				rMax = 2.0*(loop[u])[dim-2];
    			}
    			else {
    				params.a = 2.0*(loop[u])[dim-2];
    				rMax = 2.0*(loop[w])[dim-2]; 
    			}
    			rMin = params.a;
    			rGuess = (rMin + rMax)/2.0;
    			zi = 0.5*brentRootFinder(&Lfway_gsl,rGuess,rMin,rMax,tol4);
    			params.a = rL;
    			params.b = 2.0*zi;
    			ti = 0.5*( - beta + TIntegral(E,&params));
    			if (v>M/2)
    				ti *= -1.0;
    			
    			(loop[v])[dim-2] = zi;
    			(loop[v])[dim-1] = ti;
    		}
    		params.f *= 0.5;
    	}
    	
    	for (uint j=0; j<N/2; j++) {
    		loop[N/2+j] = loop[j]; 
    		loop[N/2+j] *= -1.0; // left hand side by inversion
    		//(loop[N/2+j])[dim-2] *= -1.0; // left hand side by mirroring
    	}
    
    }

	cout << "printing to " << file << endl;
	loop.save(file);
	loop.clear();

}

return 0;
}

