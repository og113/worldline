/*
	definitions for program to generate unit loops in D dimensions for use in worldline programs.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "simple.h"
#include "genloop.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - defining basic quantities
		2 - getting inputs from argv
		3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*-------------------------------------------------------------------------------------------------------------------------
	1 - defining basic quantities
		loops = 2^k
-------------------------------------------------------------------------------------------------------------------------*/

#define dim 4
uint k = 0;

/*-------------------------------------------------------------------------------------------------------------------------
	2 - getting inputs from argv
-------------------------------------------------------------------------------------------------------------------------*/
if (argc == 2) k = stringToNumber<uint>(argv[1]);
else if (argc % 2 && argc>1) {
for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("k")==0) k = stringToNumber<uint>(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

uint loops = pow(2,k);
cout << "generating " << loops << " unit loops in " << dim << " dimensions" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - initialising loops
-------------------------------------------------------------------------------------------------------------------------*/
Point<dim> p;
p[0] += 1.0;
cout << "p:" << endl << p << endl;

uint K = 2;
uint Seed = 0;
/*cout << "K: ";
cin >> K;
cout << "Seed: ";
cin >> Seed;*/

Loop<dim> l(K,Seed), m(K,Seed);
l.grow();

string file = "data/temp/loop.dat";
l.save(file);
m.load(file);

cout << "l:" << endl << l << endl;
cout << "length = " << l.checkLength() << endl;
cout << "m:" << endl << m << endl;
cout << "length = " << m.checkLength() << endl;

return 0;
}
