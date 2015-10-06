/*
	definitions for program to generate unit loops in D dimensions for use in worldline programs.
*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include "simple.h"
#include "parameters.h"
#include "genloop.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - getting parameters
		1 - initializing mpi
		2 - getting inputs from argv
		3 - starting parameter loop
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	0. getting parameters
----------------------------------------------------------------------------------------------------------------------------*/

//dimension
#define dim 4

// parameters
ParametersRange pr;
pr.load("inputs");
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	1. initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int Nw, rank, root = 0; // Nw, number of workers

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &Nw);

if (rank==root)
	cout << "starting loop with " << Nw << " nodes" << endl;

/*-------------------------------------------------------------------------------------------------------------------------
	2 - getting inputs from argv
-------------------------------------------------------------------------------------------------------------------------*/

if (argc == 2) p.K = stringToNumber<uint>(argv[1]);
else if (argc % 2 && argc>1) {
for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("k")==0 || id.compare("K")==0) p.K = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("l")==0 || id.compare("Nl")==0) p.Nl = stringToNumber<uint>(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. starting parameter loop
		- dealing with parameters
		- initializing data arrays
----------------------------------------------------------------------------------------------------------------------------*/

// parameter loops
uint Npl = 1; // number of parameter loops
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label)) {
	Npl = (pr.Steps)[label-1];
	if (rank==root)
		cout << "looping " << label << " over " << Npl << " steps" << endl;
}

// starting loop
for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pr.toStep(label) && pl>0)
		p.step(pr);
	
	uint Length = pow(2,p.K);
	uint Npw = (uint)p.Nl/Nw;

	if (rank==root) {
		// checking relevant ratios are integers
		if (p.Nl%Nw!=0) {
			cerr << "Relevant ratios are not all integers:" << endl;
			cerr << "Nl = " << p.Nl << ", Nw = " << Nw << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		cout << "generating " << p.Nl << " unit loops each of " << Length << " points in " << dim << " dimensions" << endl;
	}
		

	/*-------------------------------------------------------------------------------------------------------------------------
		4 - making and saving loops
	-------------------------------------------------------------------------------------------------------------------------*/

	string file, asciiFile;
	uint Seed = time(NULL), id;
	Loop<dim> loop(p.K,Seed);

	for (uint j=0; j<Npw; j++) {
		id = rank*Npw+j;
		file = "data/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_run_"+nts<uint>(id)+".dat";
		loop.grow();
		loop.save(file);
		loop.clear();
		loop.setSeed(time(NULL)+id*1000+3);
	}

	/*-------------------------------------------------------------------------------------------------------------------------
		4 - printing loops to graph, can only visualise output if dim=2
	-------------------------------------------------------------------------------------------------------------------------*/

	/*if (rank==root && dim==2) {
		asciiFile = "data/temp/loopAscii.dat";
		loop.grow();
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
		}
	}*/

}

MPI_Barrier(MPI_COMM_WORLD);
MPI::Finalize();

return 0;
}
