/*
	definitions for program to generate unit loops in D dimensions for use in worldline programs.
*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "simple.h"
#include "parameters.h"
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
-------------------------------------------------------------------------------------------------------------------------*/

#define dim 4

ParametersRange pr;
pr.load("inputs");
if (pr.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}
Parameters p = pr.Min;
uint Length = pow(2,p.K);

/*-------------------------------------------------------------------------------------------------------------------------
	2 - getting inputs from argv
-------------------------------------------------------------------------------------------------------------------------*/
if (argc == 2) p.K = stringToNumber<uint>(argv[1]);
else if (argc % 2 && argc>1) {
for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("k")==0 || id.compare("K")==0) p.K = stringToNumber<uint>(argv[2*j+2]);
		if (id.compare("l")==0 || id.compare("loops")==0) p.Loops = stringToNumber<uint>(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

cout << "generating " << p.Loops << " unit loops each of " << Length << " points in " << dim << " dimensions" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
string file, asciiFile;
uint Seed = time(NULL);
Loop<dim> loop(p.K,Seed);
Metropolis<dim> met(loop,p,Seed);

for (uint j=0; j<p.Loops; j++) {
	file = "data/temp/loop_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+"_run_"+nts<uint>(j)+".dat";
	loop.grow();
	if (abs(p.g)>MIN_NUMBER && p.Nms>0) {
		met.step(p.Nms*Length);
		met.setSeed(time(NULL)+j+2);
	}
	loop.save(file);
	loop.clear();
	loop.setSeed(time(NULL)+j+3);
}


//asciiFile = "data/temp/loopAscii.dat";
/*loop.grow();
//loop.saveAscii(asciiFile);
cout << "0 metropolis runs, V = " << V0(loop) << endl;
if (abs(p.g)>MIN_NUMBER) {
	for (uint j=0 ;j<8; j++) {
		uint runs = pow(10,j);
		met.setSeed(time(NULL)+j+2);
		met.step(runs);
		cout << nts<uint>(runs) << " metropolis runs, V = " << V0(loop) << endl;
		//asciiFile = "data/temp/loopAsciiMet_run_"+nts<uint>(j)+".dat";
		//loop.saveAscii(asciiFile);
	}
}*/

return 0;
}
