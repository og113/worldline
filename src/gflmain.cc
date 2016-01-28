/*
	definitions for program to generate unit fourier loops in D dimensions for use in worldline programs.
*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "simple.h"
#include "folder.h"
#include "genfloop.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - getting argv
		2 - getting parameters
		3 - checking if need to do anything
		3 - if needed, constructing s0 loops
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
	
/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

// data to print
string inputsFile = "inputs3";
bool ascii = false;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("ascii")==0) ascii = (stn<uint>(argv[2*j+2])==1);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

cout << "using inputs file " << inputsFile << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	2. getting parameters
----------------------------------------------------------------------------------------------------------------------------*/

//dimension
#define dim 4

// parameters
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3 - checking if loops already exist
-------------------------------------------------------------------------------------------------------------------------*/

bool constructLoops = false;

/*----------------------------------------------------------------------------------------------------------------------------
	4. if needed, constructing s0 loops
		- dealing with parameters
		- initializing data arrays
----------------------------------------------------------------------------------------------------------------------------*/

// parameter loops
uint Npl = 1; // number of parameter loops
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label) && label==Parameters::k) {
	Npl = (pr.Steps)[label-1];
	cout << "looping " << label << " over " << Npl << " steps" << endl;
}

// starting loop
for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pr.toStep(label) && pl>0)
		p.step(pr);
	else if (pr.toStep(label) && label==Parameters::nl)
		p.Nl = (pr.Max).Nl;

	
	// constructing folders
	FilenameAttributes faMin, faMax;
	faMin.Directory = "data/s0/floops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K);
	faMin.Timenumber = "";
	faMin.ID = "loop";
	faMin.Suffix = ".dat";
	faMax = faMin;
	(faMin.Extras).push_back(StringPair("run","0"));
	(faMax.Extras).push_back(StringPair("run",nts<uint>(p.Nl-1)));

	Folder folder(faMin,faMax);
	constructLoops = folder.size()<p.Nl;

	if (constructLoops) {

		uint Length = pow(2,p.K);
		cout << "generating " << p.Nl << " unit loops each of " << Length << " points in " << dim << " dimensions" << endl;	

		/*-------------------------------------------------------------------------------------------------------------------------
			4 - making and saving loops
		-------------------------------------------------------------------------------------------------------------------------*/

		string file;
		uint Seed = time(NULL), id;
		FLoop<dim> loop(p.K,Seed);

		for (uint j=0; j<p.Nl; j++) {
			id = j;
			file = "data/s0/floops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_run_"+nts<uint>(id)+".dat";
			loop.grow();
			loop.save(file);
			loop.clear();
			loop.setSeed(time(NULL)+id*1000+3);
		}
		
		}
	else {
		cout << "loop files already present" << endl;
	}

/*-------------------------------------------------------------------------------------------------------------------------
			4 - printing one loop in ascii
-------------------------------------------------------------------------------------------------------------------------*/

	if (ascii) {
		string asciiFile = "data/temp/asciiFloop.dat";
		uint Seed = time(NULL)+p.Nl;
		Loop<dim> loop(p.K,Seed);
		loop.grow();
		loop.saveAscii(asciiFile);
		cout << "printed loop in ascii to " << asciiFile << endl;
	}

}

return 0;
}
