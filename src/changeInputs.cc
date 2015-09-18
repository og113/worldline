/*-------------------------------------------------------------------------------------------------------------------------
 	program to change inputs or options
 -------------------------------------------------------------------------------------------------------------------------*/

#include <string>
#include <fstream>
#include <iostream>
#include "simple.h"
#include "folder.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. defining key parameters
	2. getting argv
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

/*-------------------------------------------------------------------------------------------------------------------------
	1. defining key parameters
-------------------------------------------------------------------------------------------------------------------------*/

string value;
string fi = "inputs", fo;

/*-------------------------------------------------------------------------------------------------------------------------
	2. getting argv
-------------------------------------------------------------------------------------------------------------------------*/

// getting filenames
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("f")==0 || id.compare("fi")==0 || id.compare("fileIn")==0) fi = (string)(argv[2*j+2]);
		else if (id.compare("fo")==0 || id.compare("fileOut")==0) fo = (string)(argv[2*j+2]);
	}
}
if (fo.empty()) fo = fi;

Parameters p;
p.load(fi);                   

// getting parameters
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("LoopMin")==0) p.LoopMin = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("LoopMax")==0) p.LoopMax = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("Ng")==0) p.Ng = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("Nms")==0) p.Nms = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("K")==0) p.K = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("g")==0) p.g = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("f")==0 || id.compare("fi")==0 || id.compare("fileIn")==0);
		else if (id.compare("fo")==0 || id.compare("fileOut")==0);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

// saving result
p.save(fo);

return 0;
}
