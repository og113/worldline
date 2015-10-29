/*
	binaryToAscii
		quick program to load binary files containing a vector of doubles and print ascii ones
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "simple.h"
#include "print.h"

using namespace std;

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

typedef number T;

string binaryFile = "";
string asciiFile = "";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("b")==0 || id.compare("binary")==0) binaryFile = (string)argv[2*j+2];
		else if (id.compare("a")==0 || id.compare("ascii")==0) asciiFile = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

if (binaryFile.empty() || asciiFile.empty()) {
	cerr << "must provide two files as input" << endl;
	cerr << binaryFile << ", " << asciiFile << endl;
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. loading and printing
----------------------------------------------------------------------------------------------------------------------------*/

vector<T> v;

loadVectorBinary<T>(binaryFile,v);

saveVectorAscii<T>(asciiFile,v);

cout << "saved binary contents of " << binaryFile << " to " << asciiFile << " in ascii format." << endl;

return 0;
}
