/*
	mpi test program
*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include "simple.h"
#include "folder.h"
#include "genloop.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - initializing mpi
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

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
	
cout << "node " << rank << endl;

MPI_Barrier(MPI_COMM_WORLD);
MPI::Finalize();

return 0;
}
