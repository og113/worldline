/*
	dimReduce
		program to reduce dimensions of worldlines from dim=4 to dimred=2
*/

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "analysis.h"
#include "check.h"
#include "folder.h"
#include "genloop.h"
#include "nrloop.h"
#include "print.h"
#include "simple.h"

using namespace std;

#define dim 4
#define dimred 2

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

Filename inFile, outFile;
uint K=0;

// getting argv
if (argc == 2 ) {
	inFile = (string)argv[1];
}
else {
	cerr << "must input filename of 4d loop" << endl;
	return 1;
}

// finding K
vector<number> v;
loadVectorBinary< vector<number> >(inFile,v);
K=log2(v.size()/dim); // the -1 assumes a translation invariance lagrange multiplier (removed as floor result anyway)

// assembling outFile
string dir = inFile.Directory;
size_t found = (inFile()).find((string)"dim_");
dir[found+4] = '2';
outFile = inFile;
outFile.Directory = dir;

/*----------------------------------------------------------------------------------------------------------------------------
	2. loading and printing
----------------------------------------------------------------------------------------------------------------------------*/

// defining vectors
uint N = pow(2,K);
uint zm = dim+1;
uint zmred = dimred+1;
vec x(N*dim+zm), xred(N*dimred+zmred);

// loading loop
loadVectorBinary(inFile,x);

// dimReduce
dimReduce(x, dim, N, xred, dimred, zmred);

// saving reduced loop
saveVectorBinary(outFile,xred);

cout << "dim = "<< dim << ": " << inFile << endl << "dim = " << dimred << ": " << outFile << endl;

return 0;
}
