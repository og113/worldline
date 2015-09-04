/*----------------------------------------------------------------------------------------------------------------------------
	wrapper
		program to give simple mpi wrapper for main
----------------------------------------------------------------------------------------------------------------------------*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <mpi.h>
#include "folder.h"
#include "genloop.h"
#include "evalloop.h"
#include "simple.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining loop quantitites
		2 - defining required nodes
		3 - coordinating files
		4 - evaluating loop quantitites
		5 - evaluating errors
		6 - printing results
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1. defining loop quantitites
	
	n.b. these should be put together into a parameters struct that can be saved and loaded
----------------------------------------------------------------------------------------------------------------------------*/

#define dim 2
Parameters p;
p.load("inputs");
if (p.empty()) {
	cerr << "Parameters empty, nothing in inputs file" << endl;
	return 1;
}
uint Length = pow(2,p.K);

/*----------------------------------------------------------------------------------------------------------------------------
	2. defining required nodes
		- initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int nodesCoord = 1;
int nodesWorker = 10;
int nodesTot = nodesCoord + nodesWorker;

int nodes, rank;
int returnValue = 0;

MPI::Init(argc, argv);
MPI_Comm_size(MPI_COMM_WORLD, &nodes);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if (nodes==nodesTot) {
	if (rank==0) {
		cout << "running " << nodes << " nodes" << endl;
	}
}
else {
	if (rank==0) {
		cerr << "have " << nodes << " nodes available" << endl;
		cerr << "require " << nodesTot << " nodes to run" << endl;
	}
	MPI::Finalize();
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. coordinating files
----------------------------------------------------------------------------------------------------------------------------*/

if (rank<nodesCoord) {

}

/*----------------------------------------------------------------------------------------------------------------------------
	4. evaluating loop quantitites
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	5. evaluating errors
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	6. printing results
----------------------------------------------------------------------------------------------------------------------------*/

MPI::Finalize();

return 0;
}
