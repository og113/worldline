/*
	quick test of genFloop
*/


#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>       /* time */
#include <gsl/gsl_sf_exp.h>
#include "genfloop.h"
#include "genloop.h"
#include "simple.h"


using namespace std;

int main() {

#define dim 2

uint K = 2, Seed = time(NULL)+3;
string file = "data/temp/asciiFloop.dat";

FLoop<2> fl(K,Seed);

fl.grow();

fl.saveAscii(file);
cout << "printed floop to " << file << endl;

number s0, l, errorl, v1r, errorv1r, a = 0.1, tol = 1.0e-2;
uint calls = 1e5;
s0 = S0<2>(fl);
l = L<2>(fl,errorl);
v1r = V1r<2>(fl,0.1,tol,calls,errorv1r);

cout << "a = " << a << endl;
cout << "S0 = " << s0 << endl;
cout << "L = " << l << ", error = " << errorl << endl;
cout << "V1r = " << v1r << ", error = " << errorv1r << endl;
cout << "V1r-pi*L/a = " << v1r-pi*l/a << endl;

return 0;
}
