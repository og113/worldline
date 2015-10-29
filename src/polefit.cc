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
  number * t; // variable of function
  number * y; // function to fit
  number * sigma; // sigma = 1/sqrt(weights)
};

// gsl function
int wt_f (const gsl_vector * x, void *data, gsl_vector * f) {

	size_t Nt = ((struct data *)data)->Nt;
	size_t Nk = ((struct data *)data)->Nk;
	number *t = ((struct data *)data)->t;
	number *y = ((struct data *)data)->y;
	number *sigma = ((struct data *)data)->sigma;
	number a, b, t2, Yi;

	for (uint i=0; i<Nt; i++) {
		/* Model Yi = 1 + 2T^2 Sum[a_k * b_k/(T^2-b_k^2), k ] */
		Yi = 1.0;
		for (uint j=0; j<Nk; j++) {
			a = gsl_vector_get (x,j);
			b = gsl_vector_get (x,Nk+j);
			t2 = t[i]*t[i];
			Yi += 2.0*t2*a*b/(t2-b*b);
		}
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}

	return GSL_SUCCESS;
}


// gsl function with derivatives
int wt_df (const gsl_vector * x, void *data, gsl_matrix * J) {

	size_t Nt = ((struct data *)data)->Nt;
	size_t Nk = ((struct data *)data)->Nk;
	number *t = ((struct data *)data)->t;
	number *sigma = ((struct data *)data)->sigma;
	number a, b, fract, t2;

	for (uint i=0; i<Nt; i++) {
	/* Jacobian matrix J(i,j) = dfi / dxj, */
	/* where fi = (Yi - yi)/sigma[i],      */
	/*       Yi = 1 + 2T^2 Sum[a_k * b_k/(T^2-b_k^2), k ]  */
	/* and the xj are the parameters (a_k,b_k) where k=0,...,(Nk-1) */
		for (uint j=0; j<Nk; j++) {
			a = gsl_vector_get (x,j);
			b = gsl_vector_get (x,Nk+j);
			t2 = t[i]*t[i];
			fract = 1.0/(t2-b*b);
			gsl_matrix_set (J, i, j, 	(2.0*t2*b*fract					)/sigma[i]	);
			gsl_matrix_set (J, i, Nk+j, (2.0*t2*a*fract*fract*(t2+b*b)	)/sigma[i] 	);
		}
	}
	
	return GSL_SUCCESS;
}

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

string inputsFolder = "inputs";
string dataFile = "results/s0/WvsgBT.dat";
uint Nk = 5, Nt = 0;

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
	Nt = countLines(dataFile);
	cout << "fitting data from " << dataFile << " with " << Nk << " parameters and " << Nt << " data points" << endl;
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

vector<number> T, Y, Sigma;

loadVectorAsciiColumn(dataFile,T,8);
loadVectorAsciiColumn(dataFile,Y,9);
loadVectorAsciiColumn(dataFile,Sigma,10);

if (T.size()!=Y.size() || T.size()!=Sigma.size()) {
	cerr << "polefit error: data not all same size: T.size() = " << T.size() << ", Y.size() = " << Y.size()\
				 << "Sigma.size() = " << Sigma.size() << endl;
	return 1;
}

/*-------------------------------------------------------------------------------------------------------------------------
	4 - fitting data
			to <W>(T) = 1 + Sum[a_k*T^2/(T-b_k),k]
			for k=0,...,(Nk-1)
-------------------------------------------------------------------------------------------------------------------------*/

size_t maxiter = 100, iter = 0;
number res_tol = 1.0e-8, error = 1.0;
gsl_vector *res = gsl_vector_calloc(Nt);
gsl_vector *resw = gsl_vector_calloc(Nt);
int status;

// struct for holding data
struct data d;
d.Nt = Nt;
d.Nk = Nk;
d.t = &T[0];
d.y = &Y[0];
d.sigma = &Sigma[0];

// setting fdf function
gsl_multifit_function_fdf f;
f.f = &wt_f;
f.df = &wt_df;   // set to NULL for finite-difference Jacobian 
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
gsl_vector_view sigma_gsl = gsl_vector_view_array(&Sigma[0], Nt);

// initializing solver
const gsl_multifit_fdfsolver_type * solverType = gsl_multifit_fdfsolver_lmder;
gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(solverType, Nt, 2*Nk);

// initialize starting point and weights
status = gsl_multifit_fdfsolver_set(solver, &f, &ab_gsl.vector);
if (status!=0) {
	cerr << "gsl gsl_multifit_fdfsolver_set error: " << status << endl;
	return 1;
}
res = (*solver).f;
gsl_vector_mul(resw,&sigma_gsl.vector);
error = gsl_blas_dnrm2(res);
cout << "initial residual norm: " << error << endl;

// solving for fit
/*while (error>res_tol || iter>maxiter) {
	status = gsl_multifit_fdfsolver_iterate(solver);
	if (status!=0) {
		cerr << "gsl gsl_multifit_fdfsolver_iterate error: " << status << endl;
		return 1;
	}
	res = (*solver).f;
	gsl_vector_mul(resw,&sigma_gsl.vector);
	error = gsl_blas_dnrm2(res);
	iter++;
}
if (error>res_tol) {
	cerr << "gsl fit failed to converge after " << maxiter << " iterations" << endl;
	return 1;
}
else {
	cout << "gsl fit converged after " << iter << " iterations" << endl;
	cout << "final residual norm: " << error << endl;
}
*/
// freeing solver
gsl_multifit_fdfsolver_free (solver);

/*----------------------------------------------------------------------------------------------------------------------------
	5. printing results
----------------------------------------------------------------------------------------------------------------------------*/
/*
number rate = 0.0;

for (uint k=0; k<Nk; k++) {
	rate += gsl_sf_exp(-ab[Nk+k])*ab[k]/ab[Nk+k];
}

rate *= 1.0/16.0/pi/pi;*/

return 0;
}
