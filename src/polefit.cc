/*
	polefit
		fitting to poles of <W>(T)
*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <fstream>
#include <vector>
#include "simple.h"
#include "parameters.h"
#include "print.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		A - defining functions
		
		1 - getting argv
		2 - getting parameters
		3 - loading data
		4 - fitting data
		5 - printing results
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	A. defining functions
----------------------------------------------------------------------------------------------------------------------------*/

// to hold data
struct data {
  size_t Nt;
  size_t Nk;
  number * t;
  number * w;
};

// gsl function
int wt_f (const gsl_vector * x, void *data, gsl_vector * f) {

	size_t Nt = ((struct data *)data)->Nt;
	size_t Nk = ((struct data *)data)->Nk;
	number *t = ((struct data *)data)->t;
	number *w = ((struct data *)data)->w;
	number a, b, t2, Wi;

	for (uint i=0; i<Nt; i++) {
		/* Model Wi = 1 + 2T^2 Sum[a_k * b_k/(T^2-b_k^2), k ] */
		Wi = 1.0;
		for (uint j=0; j<Nk; j++) {
			a = gsl_vector_get (x,j);
			b = gsl_vector_get (x,Nk+j);
			t2 = t[i]*t[i];
			Wi += 2.0*t2*a*b/(t2-b*b);
		}
		gsl_vector_set (f, i, Wi - w[i]);
	}

	return GSL_SUCCESS;
}


// gsl function with derivatives
int wt_df (const gsl_vector * x, void *data, gsl_matrix * J) {

	size_t Nt = ((struct data *)data)->Nt;
	size_t Nk = ((struct data *)data)->Nk;
	number *t = ((struct data *)data)->t;
	number a, b, fract, t2;

	for (uint i=0; i<Nt; i++) {
	/* Jacobian matrix J(i,j) = dfi / dxj, */
	/* where fi = (Wi - wi)/sigma[i],      */
	/*       Wi = 1 + 2T^2 Sum[a_k * b_k/(T^2-b_k^2), k ]  */
	/* and the xj are the parameters (a_k,b_k) where k=0,...,(Nk-1) */
		for (uint j=0; j<Nk; j++) {
			a = gsl_vector_get (x,j);
			b = gsl_vector_get (x,Nk+j);
			t2 = t[i]*t[i];
			fract = 1.0/(t2-b*b);
			gsl_matrix_set (J, i, j, 	2.0*t2*b*fract					);
			gsl_matrix_set (J, i, Nk+j, 2.0*t2*a*fract*fract*(t2+b*b) 	);
		}
	}
	
	return GSL_SUCCESS;
}

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

string inputsFolder = "inputs";
string dataFile = "";
uint Nk = 10, Nt = 0;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFolder = (string)argv[2*j+2];
		else if (id.compare("f")==0 || id.compare("file")==0) inputsFolder = (string)argv[2*j+2];
		else if (id.compare("Nk")==0) Nk = stn<uint>(argv[2*j+2]);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

if (dataFile.empty()) {
	cerr << "dataFile empty: ending polefit early" << endl;
	return 1;
}
else {
	cout << "fitting data from " << dataFile << " with " << Nk << " parameters" << endl;
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. getting parameters
----------------------------------------------------------------------------------------------------------------------------*/

//dimension
#define dim 4

// parameters
ParametersRange pr;
pr.load(inputsFolder);
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. loading data
----------------------------------------------------------------------------------------------------------------------------*/

string Tfile = "";
string Wfile = "";
string weightsFile = "";
vector<number> T, W, weights;

loadVectorBinary(Tfile,T);
loadVectorBinary(Wfile,W);
loadVectorBinary(weightsFile,weights);

if (T.size()==W.size() && T.size()==weights.size()) {
	Nt = T.size();
}
else {
	cerr << "polefit error: data not all same size: T.size() = " << T.size() << ", W.size() = " << W.size()\
				 << "weights.size() = " << weights.size() << endl;
	return 1;
}

/*-------------------------------------------------------------------------------------------------------------------------
	4 - fitting data
			to <W>(T) = 1 + Sum[a_k*T^2/(T-b_k),k]
			for k=0,...,(Nk-1)
-------------------------------------------------------------------------------------------------------------------------*/

size_t maxiter = 100;
number xtol = 1.0e-8;
number gtol = 1.0e-8;
number ftol = 0.0;
int status, info;

// struct for holding data
struct data d;
d.Nt = Nt;
d.Nk = Nk;
d.w = &W[0];

// setting fdf function
gsl_multifit_function_fdf f;
f.f = &wt_f;
f.df = &wt_df;   /* set to NULL for finite-difference Jacobian */
f.n = Nt;
f.p = 2*Nk;
f.params = &d;

// pick starting point
vector<number> ab_init(2*Nk,1.0);
for (uint k=0; k<Nk; k++) {
	ab_init[Nk+k] += (number)k;
}

// defining gsl vectors
gsl_vector_view ab_gsl = gsl_vector_view_array (&ab_init[0], 2*Nk);
gsl_vector_view weights_gsl = gsl_vector_view_array(&weights[0], Nt);

// initializing solver
gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmder, Nt, 2*Nk);

// initialize starting point and weights
gsl_multifit_fdfsolver_wset(solver, &f, &ab_gsl.vector, &weights_gsl.vector);

// solving for fit
status = gsl_multifit_fdfsolver_driver(solver, maxiter, xtol, gtol, ftol, &info);

// freeing solver
gsl_multifit_fdfsolver_free (solver);

/*----------------------------------------------------------------------------------------------------------------------------
	5. printing results
----------------------------------------------------------------------------------------------------------------------------*/

number rate = 0.0;

for (uint k=0; k<Nk; k++) {
	rate += gsl_sf_exp(-ab[Nk+k])*ab[k]/ab[Nk+k];
}

rate *= 1.0/16.0/pi/pi;

return 0;
}
