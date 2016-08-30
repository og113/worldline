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

// getting parameters
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("Nl")==0) p.Nl = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("Ng")==0) p.Ng = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("Nig")==0) p.Nig = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("Nsw")==0) p.Nsw = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("Npsw")==0) p.Npsw = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("K")==0) p.K = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("G")==0) p.G = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("B")==0) p.B = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("T")==0) p.T = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("Epsi")==0) p.Epsi = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("Mu")==0) p.Mu = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("P1")==0) p.P1 = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("P2")==0) p.P2 = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("P3")==0) p.P3 = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("P4")==0) p.P4 = stringToNumber<number>(argv[2*j+2]);
		else if (id.compare("Lambda")==0) p.Lambda = stringToNumber<number>(argv[2*j+2]);
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
