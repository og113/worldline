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

#define dim 4

uint K = 4, Seed = time(NULL)+3;
string afile = "data/temp/asciiFloop.dat", bfile = "data/temp/binaryFloop.dat", fcfile = "temp/fc.dat";

FLoop<dim> fl(K,Seed);

fl.grow();

ofstream os;
os.open(fcfile.c_str());
FCoeff<dim> fc;
fc.writeBinary(os);
os.close();
cout << "size of FCoeff = " << countDoubles(fcfile) << endl;

fl.save(bfile);

FLoop<dim> fl2(K,Seed);
fl2.load(bfile);

os.open(fcfile.c_str());
(fl[0]).writeBinary(os);
os.close();
cout << "size of FCoeff = " << countDoubles(fcfile) << endl;

fl.saveAscii(afile);
cout << "printed floop to " << afile << endl;

number s0, l, errorl, v1r, errorv1r, a = 0.1, tol = 1.0e-2;
uint calls = 1e5;
s0 = S0<dim>(fl);
l = L<dim>(fl,errorl);
v1r = V1r<dim>(fl,0.1,tol,calls,errorv1r);

cout << "a = " << a << endl;
cout << "S0 = " << s0 << endl;
cout << "L = " << l << ", error = " << errorl << endl;
cout << "V1r = " << v1r << ", error = " << errorv1r << endl;
cout << "V1r-pi*L/a = " << v1r-pi*l/a << endl;

return 0;
}
