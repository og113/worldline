/*----------------------------------------------------------------------------------------------------------------------------
	nrmpi
		mpi wrapper for nrmain
----------------------------------------------------------------------------------------------------------------------------*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include "simple.h"
#include "parameters.h"
#include "nrmain_fn.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - initializing mpi
		1 - sorting inputs
		2 - running nrmain
		3 - clearing up
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

using namespace std;

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	0. initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int Nw, rank, root = 0; // Nw, number of workers

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &Nw);

if (rank==root)
	cout << "starting glmain with " << Nw << " nodes" << endl;	

/*----------------------------------------------------------------------------------------------------------------------------
	1. sorting inputs
----------------------------------------------------------------------------------------------------------------------------*/

// argc and argv
uint argcSerial = 1;
string mpiInputsFile = "nrinputs/mpi/inputs0";
vector<string> argvSerial(1);
argvSerial[0] = "./nrmpi";

// printing argv
if (rank==root) {
	cout << "input argv to ./nrmpi:" << endl;
	for (int j=0; j<argc; j++) {
		cout << argv[j];
		if (j<(argc-1))
			cout << " ";
		else
			cout << endl;
	}
}

// getting mpi and serial argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("mpiinputs")==0 || id.compare("mpiInputs")==0 || id.compare("inputs")==0)
			mpiInputsFile = (string)argv[2*j+2];
		else {
			argvSerial.push_back((string)argv[2*j+1]);
			argvSerial.push_back((string)argv[2*j+2]);
			argcSerial += 2;
		}
	}
}
else if (rank==root) {
	cerr << "must provide an even number of arguments after ./nrmpi" << endl;
	MPI_Abort(MPI_COMM_WORLD,1);
}
	
sleep(rank*2);
	
// loading parameters
ParametersRange prMpi;
prMpi.load(mpiInputsFile);
if (rank==root) {
	Parameters pMpi = prMpi.Min;
	if (pMpi.empty()) {
		cerr << "Parameters empty: nothing in inputs file: " << mpiInputsFile << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if (Nw>(int)prMpi.totalSteps()) {
		cerr << "Error: too few steps, " << prMpi.totalSteps() << ", for nodes, " << Nw << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
}

// distributing parameters
string inputsSerial = mpiInputsFile+"_rank_"+nts(rank);
ParametersRange prSerial = prMpi;
prSerial.Min = prMpi.position(rank);
uint j=0;
bool foundToStep = false;
while (j<Parameters::Size && !foundToStep) {
	if (prMpi.toStep((Parameters::Label)(j+1))) {
		(prSerial.Steps)[j] = 0;
		if (rank==root && (int)(prMpi.Steps)[j]!=Nw) {
			cerr << "Error: number of steps on first parameter, " << (prMpi.Steps)[j];
			cerr << ", not equal to number of nodes, " << Nw << endl;
		}
		foundToStep = true;
	}
	j++;
}
prSerial.save(inputsSerial);

// adding inputs filename to argv
argvSerial.push_back("-inputs");
argvSerial.push_back(inputsSerial);
argcSerial += 2;

/*----------------------------------------------------------------------------------------------------------------------------
	2. running nrmain
----------------------------------------------------------------------------------------------------------------------------*/

uint returnValue = 0;
returnValue = nrmain_fn(argcSerial,argvSerial);

if (returnValue!=0) {
	cerr << "----------------------------------------------------------------------" << endl;
	cerr << "return " << returnValue << " for node " << rank << " on running nrmain" << endl;
	cerr << "----------------------------------------------------------------------" << endl;
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. clearing up
----------------------------------------------------------------------------------------------------------------------------*/

MPI::Finalize();

return 0;
}
