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

#define dim 2

Parameters p;
/*p.LoopMin = 1;
p.LoopMax = 1e4;
p.Ng = 1e3;
p.K = 5;
p.g = 0.0;*/
p.load("inputs");
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}
uint Loops = p.LoopMax-p.LoopMin+1; // number of loops
uint Length = pow(2,p.K);

/*-------------------------------------------------------------------------------------------------------------------------
	2 - getting inputs from argv
-------------------------------------------------------------------------------------------------------------------------*/
if (argc == 2) p.K = stringToNumber<uint>(argv[1]);
else if (argc % 2 && argc>1) {
for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("k")==0 || id.compare("K")==0) p.K = stringToNumber<uint>(argv[2*j+2]);
		if (id.compare("l")==0 || id.compare("loops")==0) p.K = stringToNumber<uint>(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

cout << "generating " << Loops << " unit loops each of " << Length << " points in " << dim << " dimensions" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	3 - making and saving loops
-------------------------------------------------------------------------------------------------------------------------*/
string file, asciiFile;
uint Seed = time(NULL);
Loop<dim> loop(p.K,Seed);

for (uint j=p.LoopMin; j<=p.LoopMax; j++) {
file = "data/temp/loop_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+"_run_"+nts<uint>(j)+".dat";
loop.grow();
loop.save(file);
/*if (abs(loop.checkLength()-1.0)>MIN_NUMBER*Length)
	cerr << "loop error: length = " << loop.checkLength() << endl;*/
loop.clear();
//Seed = time(NULL)+j;
//loop.setSeed(Seed);
}

asciiFile = "data/temp/loopAscii.dat";
Seed += 2.0;
loop.grow();
loop.saveAscii(asciiFile);
Metropolis<dim> met(loop,Seed);
met.step(1000);
asciiFile = "data/temp/loopAsciiMet.dat";
loop.saveAscii(asciiFile);

return 0;
}
