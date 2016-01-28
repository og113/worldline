/*
	circle
		program to generate circlular loop
*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf_trig.h>
#include <fstream>
#include "simple.h"
#include "parameters.h"
#include "genloop.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - getting inputs from argv
		2 - defining basic quantities
		3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

// data to print
string inputsFile = "inputs";
bool extend = false;
bool lemon = false;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("extend")==0) extend = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("lemon")==0) lemon = (stn<uint>(argv[2*j+2])!=0);
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
if (pr.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}
Parameters p = pr.Min;

// parameter loops
uint Npl = 1; // number of parameter loops
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label)) {
	Npl = (pr.Steps)[label-1];
	cout << "looping " << label << " over " << Npl << " steps" << endl;
}

// starting loop
for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pr.toStep(label) && pl>0)
		p.step(pr);
	else if (pr.toStep(label) && label==Parameters::nl)
		p.Nl = (pr.Max).Nl;
		
	uint N = pow(2,p.K);

	cout << "generating " << p.Nl << " ";
	if (!lemon)
		cout << "circular";
	else
		cout << "lemon shaped";
	cout << " unit loops each of " << N << " points in " << dim << " dimensions" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
	string file, asciiFile;
	uint Seed = time(NULL);
	Loop<dim> loop(p.K,Seed);
	Metropolis<dim> met(loop,p,Seed);
	Point<dim> p0, point;
	number R = 1.0, R0 = 1.0, M = abs(p.P4);
	if (abs(p.G)>MIN_NUMBER && abs(p.B)>MIN_NUMBER) {
		R =  1.0/p.G/p.B;
		R0 = R;
	}
	if (abs(p.Epsi)>MIN_NUMBER && extend)
		R += p.Epsi;
	
	for (uint j=0; j<p.Nl; j++) {
		if (!lemon)
			file = "data/circle/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R0)\
															+"_rank_"+nts(j)+".dat";
		else
			file = "data/lemon/loops/dim_"+nts<uint>(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R0)\
										+"_M_"+nts(M)+"_rank_"+nts<uint>(j)+".dat";
		if (j==0)
			cout << "printing to " << file << ", with runs from 0..." << p.Nl << endl;
			
		if (!lemon)	{
			number w = 2.0*PI/(number)N;
			for (uint k=0; k<N; k++) {
				point = p0;
				point[2] += R*gsl_sf_cos(w*k);
				point[3] += R*gsl_sf_sin(w*k);
				loop[k] = point;
			}
		}
		else {
			number x0 = R*p.P4/2.0;
			number rho = R*sqrt(1.0-p.P4*p.P4/4.0), angle;
			Point<dim> pl, pr;
			pl[2] -= x0;
			pr[2] += x0;
			number w_max = asin(rho/R);
			
			for (uint k=0; k<N; k++) {
				if (k<N/2) {
					angle = -w_max + (4.0*w_max/(number)N)*k;
					point = pl;
					point[2] += R*gsl_sf_cos(angle);
					point[3] += R*gsl_sf_sin(angle);
					loop[k] = point;
				}
				else {
					angle = -w_max + (4.0*w_max/(number)N)*(k-(number)N/2.0);
					point = pr;
					point[2] += -R*gsl_sf_cos(angle);
					point[3] += -R*gsl_sf_sin(angle);
					loop[k] = point;
				}
			}
		}
		
		loop.save(file);
		loop.clear();
	}

}

return 0;
}
