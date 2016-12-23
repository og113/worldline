/*
	circle
		program to generate circlular loop
*/

#include <ctime>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf_trig.h>
#include <fstream>
#include "simple.h"
#include "parameters.h"
#include "genloop.h"
#include "folder.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - getting inputs from argv
		2 - defining basic quantities
		3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

struct ShapeOptions {
	enum Option { circle, lemon, straightDisjoint, cosDisjoint};
};

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

// data to print
string inputsFile = "inputs4";
bool extend = false;
bool higherOrder = false;
bool circle = false;
bool lemon = false;
bool straightDisjoint = false;
bool cosDisjoint = false;
string shape = "circle";
string baseFolder = "";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("base")==0) baseFolder = (string)argv[2*j+2];
		else if (id.compare("extend")==0) extend = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("lemon")==0) lemon = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("circle")==0) circle = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("higherOrder")==0) higherOrder = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("straightDisjoint")==0) straightDisjoint = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("cosDisjoint")==0) cosDisjoint = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("shape")==0) shape = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

ShapeOptions::Option so = ShapeOptions::circle;
if (!shape.empty() || circle || lemon || straightDisjoint || cosDisjoint) {
	if (shape.compare("circle")==0 || circle) {
		so =ShapeOptions::circle;
		shape = "circle";
	}
	else if (shape.compare("lemon")==0 || lemon) {
		so = ShapeOptions::lemon;
		shape = "lemon";
	}
	else if (shape.compare("straightDisjoint")==0 || straightDisjoint) {
		so = ShapeOptions::straightDisjoint;
		shape = "straightDisjoint";
	}
	else if (shape.compare("cosDisjoint")==0 || cosDisjoint) {
		so = ShapeOptions::cosDisjoint;
		shape = "cosDisjoint";
	}
	else {
		cerr << "shape options not understood: " << shape << endl;
		return 1;
	}
}

cout << "using inputs file " << inputsFile << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	2 - defining basic quantities
-------------------------------------------------------------------------------------------------------------------------*/

#define dim 2

ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}

// parameter loops
uint Npl = pr.totalSteps();
cout << "looping over " << Npl << " steps" << endl;


// starting loop
for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pl>0)
		p = pr.position(pl);
		
	uint N = pow(2,p.K);

	cout << "generating " << p.Nl << " " << shape << " shaped unit loops each of ";
	cout << N << " points in " << dim << " dimensions" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
	Filename file;
	uint Seed = time(NULL);
	Loop<dim> loop(p.K,Seed);
	Metropolis<dim> met(loop,p,Seed);
	Point<dim> p0, point;
	number R = 1.0, R0 = 1.0;
	number beta = ((p.T)>sqrt(MIN_NUMBER)? 1.0/(p.T): 1.0/sqrt(MIN_NUMBER));
	number kappa = pow(p.G,3)*p.B;
	// scaling loops by g*B, so R = 1.0 always
	/*if (abs(p.G)>MIN_NUMBER && abs(p.B)>MIN_NUMBER) {
		R =  1.0/p.G/p.B;
		R0 = R;
	}*/
	if (abs(p.Epsi)>MIN_NUMBER && extend)
		R += p.Epsi;
	
	
	for (uint j=0; j<p.Nl; j++) {

		if (so==ShapeOptions::circle)	{
			file = baseFolder+"data/"+shape+"/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R0)\
															+"_rank_"+nts(j)+".dat";
			number w = 2.0*PI/(number)N, angle;
			for (uint k=0; k<N; k++) {
				point = p0;
				angle = w*k - PI/2.0;
				point[dim-2] += R*gsl_sf_cos(angle);
				point[dim-1] += R*gsl_sf_sin(angle);
				loop[k] = point;
			}
		}
		else if (so==ShapeOptions::lemon) {
			file = baseFolder+"data/"+shape+"/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R0)\
									+"_E_"+nts(p.P4)+"_rank_"+nts<uint>(j)+".dat";
			number x0 = R*p.P4/2.0;
			number rho = R*sqrt(1.0-p.P4*p.P4/4.0), angle;
			Point<dim> pl, pr;
			pl[dim-2] -= x0;
			pr[dim-2] += x0;
			number w_max = asin(rho/R);
			
			for (uint k=0; k<N; k++) {
				if (k<N/2) {
					angle = -w_max + (4.0*w_max/(number)N)*k;
					point = pl;
					point[dim-2] += R*cos(angle);
					point[dim-1] += R*sin(angle);
					loop[k] = point;
				}
				else {
					angle = -w_max + (4.0*w_max/(number)N)*(k-(number)N/2.0);
					point = pr;
					point[dim-2] += -R*cos(angle);
					point[dim-1] += -R*sin(angle);
					loop[k] = point;
				}
			}
		}
		else if (so==ShapeOptions::straightDisjoint) {
			file = baseFolder+"data/"+shape+"/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(kappa)\
					+"_T_"+nts(p.T)+"_rank_"+nts<uint>(j)+".dat";
			if (extend)	
				(file.Extras).push_back(StringPair("Lambda",nts(p.Lambda)));
			number r = sqrt(kappa/4.0/PI);
			if (higherOrder)
				r += - 3.0*sqrt(PI/4.0/kappa)*pow(p.Epsi,2)\
						- 15.0*pow(PI/kappa,3.0/2.0)*pow(p.Epsi,4)/4.0;
			if (extend)
				r *= (1.0 + p.Lambda);
			number dt = 2.0*beta/(number)N;
			for (uint k=0; k<N; k++) {
				point = p0;
				if (k<N/2) {
					point[dim-2] += r/2.0;
					point[dim-1] += -beta/2.0 + dt/2.0 + dt*k;
				}
				else {
					point[dim-2] += -r/2.0;
					point[dim-1] += beta/2.0 - dt/2.0 - dt*(k-N/2);
				}
				loop[k] = point;
			}
		}
		else if (so==ShapeOptions::cosDisjoint) {
			file = baseFolder+"data/"+shape+"/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(kappa)\
					+"_T_"+nts(p.T)+"_mu_"+nts(p.Mu)+"_lambda_"+nts(p.Lambda)+"_rank_"+nts(j)+".dat";	
			number r = sqrt(kappa/4.0/PI);
			if (higherOrder)
				r += - 3.0*sqrt(PI/4.0/kappa)*pow(p.Epsi,2)\
						- 15.0*pow(PI/kappa,3.0/2.0)*pow(p.Epsi,4)/4.0;
			number dt = 2.0*beta/(number)N;
			number w = 4.0*PI/(number)N;
			for (uint k=0; k<N; k++) {
				point = p0;
				if (k<N/2) {
					point[dim-2] += (r/2.0)*(1.0 + p.Mu*cos(w*(k+0.5)) + p.Lambda*cos(2.0*w*(k+0.5)));
					point[dim-1] += -beta/2.0 + dt/2.0 + dt*k;
				}
				else {
					point[dim-2] += -(r/2.0)*(1.0 + p.Mu*cos(w*(k+0.5)) + p.Lambda*cos(2.0*w*(k+0.5)));
					point[dim-1] += beta/2.0 - dt/2.0 - dt*(k-N/2);
				}
				loop[k] = point;
			}
		}
		
		if (j==0)
			cout << "printing to " << file << ", with runs from 0..." << p.Nl << endl;
		
		loop.save(file);
		loop.clear();
	}

}

return 0;
}
