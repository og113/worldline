/*----------------------------------------------------------------------------------------------------------------------------
	wrapper
		program to give simple mpi wrapper for main
----------------------------------------------------------------------------------------------------------------------------*/

#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <mpi.h>
#include <gsl/gsl_sf_exp.h>
#include "folder.h"
#include "genloop.h"
#include "parameters.h"
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
	
----------------------------------------------------------------------------------------------------------------------------*/

#define dim 4
Parameters p;
p.load("inputs");
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}
uint Length = pow(2,p.K);
uint dof = dim*(Length-1);

vector<number> dataS0(p.Ng,0.0), dataS02(p.Ng,0.0); // dataZ(p.Ng,0.0)
number aprxS0 = 0.0, aprxS02 = 0.0, error = 0.0; // aprxZ = 0.0

/*----------------------------------------------------------------------------------------------------------------------------
	2. defining required nodes
		- checking all the ratios are integers
		- initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int Nw = 1;
int nodesTot = 1 + Nw;
uint Nl, Npg, Npw;
//number T = 0.25;

if ((p.LoopMax-p.LoopMin)<=0) {
	cerr << "Parameters error: LoopMax<=LoopMin" << endl;
	return 1;
}
else {
	Nl = p.LoopMax-p.LoopMin+1;
}

if (Nl%p.Ng!=0 || Nl%Nw!=0 || p.Ng%Nw!=0) {
	cerr << "Ratios are not all integers:" << endl;
	cerr << "Nl = " << Nl << ", Ng = " << p.Ng << ", Nw = " << Nw << endl;
	return 1;
}
else {
	Npg = (uint)Nl/p.Ng;
	Npw = (uint)Nl/Nw;
}

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

uint loopMin, loopMax;

if (rank==0) {
	for (int k=0; k<Nw; k++) {
		loopMin = k*Npw;
		loopMax = (k+1)*Npw-1;
		
		MPI::COMM_WORLD.Send(&loopMin, 1, MPI::UNSIGNED, k+1, 0);		
		MPI::COMM_WORLD.Send(&loopMax, 1, MPI::UNSIGNED, k+1, 1);
		
		//cout << "process " << 0 << " sent " << loopMin << " and " << loopMax << " to process " << k+1 << endl;
		
	}
}
else {
	MPI::Status status;
	MPI::COMM_WORLD.Probe(0, 0, status);
	MPI::COMM_WORLD.Recv(&loopMin, 1, MPI::UNSIGNED, 0, 0, status);	
	MPI::COMM_WORLD.Probe(0, 1, status);
	MPI::COMM_WORLD.Recv(&loopMax, 1, MPI::UNSIGNED, 0, 1, status);
	
	//cout << "process " << rank << " recieved " << loopMin << " and " << loopMax << " from process 0" << endl;
}

/*----------------------------------------------------------------------------------------------------------------------------
	4. evaluating loop quantitites
----------------------------------------------------------------------------------------------------------------------------*/

if (rank>0) {
	FilenameAttributes faMin, faMax;

	faMin.Directory = "data/temp";
	faMin.Timenumber = "";
	(faMin.Extras).push_back(StringPair("dim",nts<uint>(dim)));
	(faMin.Extras).push_back(StringPair("K",nts<uint>(p.K)));
	faMax = faMin;
	(faMin.Extras).push_back(StringPair("run",nts<uint>(loopMin)));
	(faMax.Extras).push_back(StringPair("run",nts<uint>(loopMax)));

	Folder folder(faMin,faMax);
	if (folder.size()!=(loopMax-loopMin+1)) {
		cerr << "error for processor " << rank << ":" << endl;
		cerr << "folder.size() = " << folder.size() << ", loopMax-loopMin+1 = " << loopMax-loopMin+1 << endl;
		return 1;
	}
	
	uint Seed = time(NULL)+rank+1;
	Loop<dim> l(p.K,Seed);
	uint counter = 0;
	uint id;
	number s0;
	//number e_s0;
	number sums[2];

	for (uint j=loopMin; j<=loopMax; j++) {
		counter++;
		l.load(folder[j]);
	
		s0 = S0(l);
		//e_s0 = gsl_sf_exp(-s0);
		sums[0] += s0;
		sums[1] += s0*s0;
		//sums[2] += s0*s0*e_s0;
		
		if (counter==Npg) {
			id = (rank-1)*(p.Ng/Nw) + ((j+1)/Npg-1);		
			MPI::COMM_WORLD.Send(&sums, 2, MPI::DOUBLE, 0, id);
			//cout << "process " << rank << " sent message " << id << " to " << 0 << endl;
			memset(sums,0,sizeof(sums));
			counter = 0;
		}
	}
}
else { // rank==0
	number buf[2];
	MPI::Status status;
	uint count=0;
	while (count<p.Ng) {
		MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, status);
		MPI::COMM_WORLD.Recv(buf, 3, MPI::DOUBLE, status.Get_source(), status.Get_tag(), status);
		dataS0[status.Get_tag()] = buf[0];
		dataS02[status.Get_tag()] = buf[1];
		//dataZ[status.Get_tag()] = buf[2];
		//cout << "process " << rank << " recieved message " << status.Get_tag() << " from " << status.Get_source() << endl;
		count++;
		aprxS0 += buf[0];
		aprxS02 += buf[1];
		//aprxZ += buf[2];
	}
	aprxS0 /= (double)Nl;
	aprxS02 /= (double)Nl;
	// better to do immediate probing 'iprobe' so that the 0 processor can do other tasks while waiting, like computing errors
}

/*----------------------------------------------------------------------------------------------------------------------------
	5. evaluating errors
----------------------------------------------------------------------------------------------------------------------------*/

if (rank==0) {
	number aprxS0_g[p.Ng];
	number denom = (number)p.Ng*((number)p.Ng-1.0);
	for (uint j=0; j<p.Ng; j++) {
		aprxS0_g[j] = dataS0[j]/(number)Npg;
		error += pow(aprxS0_g[j]-aprxS0,2.0)/denom;
	}
	error = sqrt(error);
	
	number analytic = 0.5*(number)dof;
	number absError = (analytic-aprxS0)/analytic;
	number variance = aprxS02-aprxS0*aprxS0;
	number abs2Error = (variance-analytic)/analytic;

/*----------------------------------------------------------------------------------------------------------------------------
	6. printing results
----------------------------------------------------------------------------------------------------------------------------*/

	string timenumber = currentDateTime();	
	
	Filename rf = "results/"+timenumber+"loopGroups_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+".dat";
	FILE * ros;
	ros = fopen(((string)rf).c_str(),"w");
	for (uint j=0; j<p.Ng; j++) {
		fprintf(ros,"%12s%5i%5i%5i%5i%5i%13.5g%13.5g%13.5g\n",\
				timenumber.c_str(),dim,p.K,Nl,p.Ng,j,aprxS0_g[j],aprxS0,error);
	}
	fclose(ros);	
	cout << "results printed to " << rf << endl;
	
	rf.Timenumber = "";
	rf.ID = "loop";
	ros = fopen(((string)rf).c_str(),"a");
		fprintf(ros,"%12s%5i%5i%5i%5i%13.5g%13.5g\n",\
				timenumber.c_str(),dim,p.K,Nl,p.Ng,aprxS0,error);
	fclose(ros);	
	cout << "results printed to " << rf << endl;

	cout << "timenumber: " << timenumber << endl;
	printf("\n");
	printf("%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s\n","dim","Nl","Ng","K","S0","S02","var",\
		"error","absError","abs2Error");
	printf("%8i%8i%8i%8i%12.3g%12.3g%12.3g%12.3g%12.3g%12.3g\n",\
		dim,Nl,p.Ng,p.K,aprxS0,aprxS02,variance,error,absError,abs2Error);
	printf("\n");
}

MPI::Finalize();

return 0;
}
