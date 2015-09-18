/*----------------------------------------------------------------------------------------------------------------------------
	wrapper
		program to give simple mpi wrapper for main
----------------------------------------------------------------------------------------------------------------------------*/

#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <mpi.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
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

vector<number> dataS0(p.Ng,0.0), dataS02(p.Ng,0.0), dataV(p.Ng,0.0), dataV2(p.Ng,0.0);
number aprxS0 = 0.0, aprxS02 = 0.0, errorS0 = 0.0, aprxV = 0.0, aprxV2 = 0.0, errorV = 0.0;

/*----------------------------------------------------------------------------------------------------------------------------
	2. defining required nodes
		- checking all the ratios are integers
		- initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int Nw = 5;
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
	
	uint Seed = time(NULL)+rank+2;
	Loop<dim> l(p.K,Seed);
	uint counter = 0;
	uint id;
	number s0, v, w;
	number sums[4];

	for (uint j=0; j<Npw; j++) {
		counter++;
		l.load(folder[j]);
	
		s0 = S0(l);
		v = V0(l);
		w = gsl_sf_cos(p.g*I0(l));
		sums[0] += s0;
		sums[1] += s0*s0;
		sums[2] += w;
		sums[3] += w*w;
		
		if (counter==Npg) {
			id = (rank-1)*(p.Ng/Nw) + ((j+1)/Npg-1);		
			MPI::COMM_WORLD.Send(&sums, 4, MPI::DOUBLE, 0, id);
			//cout << "process " << rank << " sent message " << id << " to " << 0 << endl;
			memset(sums,0,sizeof(sums));
			counter = 0;
		}
	}
}
else { // rank==0
	number buf[4];
	MPI::Status status;
	uint count=0;
	while (count<p.Ng) {
		MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, status);
		MPI::COMM_WORLD.Recv(buf, 4, MPI::DOUBLE, status.Get_source(), status.Get_tag(), status);
		dataS0[status.Get_tag()] = buf[0];
		dataS02[status.Get_tag()] = buf[1];
		dataV[status.Get_tag()] = buf[2];
		dataV2[status.Get_tag()] = buf[3];
		//cout << "process " << rank << " recieved message " << status.Get_tag() << " from " << status.Get_source() << endl;
		count++;
		aprxS0 += buf[0];
		aprxS02 += buf[1];
		aprxV += buf[2];
		aprxV2 += buf[3];
	}
	aprxS0 /= (number)Nl;
	aprxS02 /= (number)Nl;
	aprxV /= (number)Nl;
	aprxV2 /= (number)Nl;
	// better to do immediate probing 'iprobe' so that the 0 processor can do other tasks while waiting, like computing errors
}

/*----------------------------------------------------------------------------------------------------------------------------
	5. evaluating errors
----------------------------------------------------------------------------------------------------------------------------*/

if (rank==0) {
	number aprxS0_g[p.Ng], aprxV_g[p.Ng];
	number denom = (number)p.Ng*((number)p.Ng-1.0);
	for (uint j=0; j<p.Ng; j++) {
		aprxS0_g[j] = dataS0[j]/(number)Npg;
		aprxV_g[j] = dataV[j]/(number)Npg;
		errorS0 += pow(aprxS0_g[j]-aprxS0,2.0)/denom;
		errorV += pow(aprxV_g[j]-aprxV,2.0)/denom;
	}
	errorS0 = sqrt(errorS0);
	errorV = sqrt(errorV);
	
	number varianceS0 = aprxS02-aprxS0*aprxS0;
	number varianceV = aprxV2-aprxV*aprxV;

/*----------------------------------------------------------------------------------------------------------------------------
	6. printing results
----------------------------------------------------------------------------------------------------------------------------*/

	string timenumber = currentDateTime();	
	
	Filename rf = "results/"+timenumber+"loopGroups_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+".dat";
	FILE * ros;
	ros = fopen(((string)rf).c_str(),"w");
	for (uint j=0; j<p.Ng; j++) {
		fprintf(ros,"%12s%5i%5i%8i%8i%10.2g%8i%13.5g%13.5g%13.5g%13.5g%13.5g%13.5g\n",\
				timenumber.c_str(),dim,p.K,Nl,p.Ng,p.g,j,aprxS0_g[j],aprxS0,errorS0,aprxV_g[j],aprxV,errorV);
	}
	fclose(ros);	
	cout << "results printed to " << rf << endl;
	
	rf = "results/loop_dim_"+nts<uint>(dim)+".dat";
	ros = fopen(((string)rf).c_str(),"a");
		fprintf(ros,"%12s%5i%5i%8i%8i%10.2g%13.5g%13.5g%13.5g%13.5g\n",\
				timenumber.c_str(),dim,p.K,Nl,p.Ng,p.g,aprxS0,errorS0,aprxV,errorV);
	fclose(ros);	
	cout << "results printed to " << rf << endl;

	cout << "timenumber: " << timenumber << endl;
	printf("\n");
	printf("%8s%8s%8s%8s%10s%12s%12s%12s%12s%12s%12s%12s%12s\n","dim","Nl","Ng","K","g","S0","S02","varS0",\
		"errorS0","V","V2","varV","errorV");
	printf("%8i%8i%8i%8i%10.2g%12.3g%12.3g%12.3g%12.3g%12.3g%12.3g%12.3g%12.3g\n",\
		dim,Nl,p.Ng,p.K,p.g,aprxS0,aprxS02,varianceS0,errorS0,aprxV,aprxV2,varianceV,errorV);
	printf("\n");
}

MPI::Finalize();

return 0;
}
