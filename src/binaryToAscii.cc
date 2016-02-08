/*
	binaryToAscii
		quick program to load binary files containing a vector of doubles and print ascii ones
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "simple.h"
#include "folder.h"
#include "genloop.h"
#include "print.h"

using namespace std;

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

Filename binaryFile, asciiFile;
bool loop = false;
uint K=0;

//dimension
#define dim 4

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("b")==0 || id.compare("binary")==0 || id.compare("fi")==0) binaryFile = (string)argv[2*j+2];
		else if (id.compare("a")==0 || id.compare("ascii")==0 || id.compare("fo")==0) asciiFile = (string)argv[2*j+2];
		else if (id.compare("l")==0 || id.compare("loop")==0) loop = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("K")==0) K = stn<uint>(argv[2*j+2]);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

if (loop && K==0) {
	cout << "input K of loop: ";
	cin >> K;
}

if (binaryFile.empty() || asciiFile.empty()) {
	cerr << "must provide two files as input" << endl;
	cerr << binaryFile << ", " << asciiFile << endl;
	return 1;
}

if (!binaryFile.exists()) {
	cerr << "binary file " << binaryFile << " doesn't exist" << endl;
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. loading and printing
----------------------------------------------------------------------------------------------------------------------------*/

vector<number> v;

if (loop) {
	Loop<dim> l(K,0);
	l.load(binaryFile);
	l.saveAscii(asciiFile);
}
else {
	loadVectorBinary< vector<number> >(binaryFile,v);
	if (v.size()>pow(2,K) && K!=0)
		v.resize(pow(2,K));
	saveVectorAscii< vector<number> >(asciiFile,v);
}

cout << "saved binary contents of " << binaryFile << " to " << asciiFile << " in ascii format." << endl;

return 0;
}
