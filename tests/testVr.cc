/*
	test of Vr and potentially other functions of loops
*/

#include <cmath>
#include <fstream>
#include <gsl/gsl_sf_exp.h>
#include <iostream>
#include <string>
#include <vector>
#include "genloop.h"
#include "simple.h"
#include "print.h"

using namespace std;

/*------------------------------------------------------------------------------------------------------------------------*/
int main() {

cout << "testVr" << endl << endl;

// dimension
#define dim 4

// parameters
Parameters p;
p.K = 8;
p.G = 1.0;
p.Epsi = 0.001;

// file
string lf = "data/s0/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_run_0.dat";

// initializing loop
uint Seed = 2;
Loop<dim> loop(p.K,Seed);
loop.load(lf);

cout << "loop file:" << lf << endl;
cout << "test of potentials:" << endl;

cout << "V0:  " << V0(loop) << endl;
cout << "V1:  " << V1(loop) << endl;
cout << "V1r: " << V1r(loop,p.Epsi) << endl;
cout << "V2r: " << V2r(loop,p.Epsi) << endl << endl;

// looPIng over a
uint Loops = 16;
number a;
cout << "looPIng over V1r(a): " << endl;
for (uint j=0; j<Loops; j++) {
	a = gsl_sf_exp(-(j+1.0));
	cout << "V1r(" << a << "): " << V1r(loop,a) << endl;
}

// looPIng over N
Loops = 8;
uint K = 4;
cout << endl << "looPIng over V1r(K): " << endl;
for (uint j=0; j<Loops; j++) {
	Loop<dim> loopK(K,Seed);
	lf = "data/s0/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(K)+"/loop_run_0.dat";
	loopK.load(lf);
	cout << "V1r(" << K << "): " << V1r(loopK,p.Epsi) << endl;
	K++;
}


return 0;
}
