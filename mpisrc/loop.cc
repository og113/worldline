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
		3 - starting parameter loop
		4 - coordinating files
		5 - evaluating loop quantitites
		6 - evaluating errors
		7 - printing results
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1. defining basic quantitites
	
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

// parameter loops
uint Npl = 1;
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label)) {
	Npl = (pr.Steps)[label];
}

// quantities to calculate
// elements are in pairs of (X,X^2): S0, V, W
// separate vector for error
uint Nq = 3;
vector<number> sums(2*Nq,0.0);
vector<number> error(Nq,0.0);

// distribution of quantity
//vector<number> dataS0(p.Ng,0.0);
//vector<number> dataS02(p.Ng,0.0);

/*----------------------------------------------------------------------------------------------------------------------------
	2. defining required nodes
		- initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int Nw, rank, root=0;
uint Nl, Npg, Npw;

MPI_Init(NULL, NULL);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &Nw);

/*----------------------------------------------------------------------------------------------------------------------------
	3. starting parameter loop
----------------------------------------------------------------------------------------------------------------------------*/

for (uint pl=0; pl<Npl; pl++) {
	if (rank==root) {
		// stepping parameters
		if (pr.toStep(label))
			p.step(pr);
		Nl = p.Loops;
		
		// checking relevant ratios are integers
		if (Nl%p.Ng!=0 || Nl%Nw!=0 || p.Ng%Nw!=0) {
			cerr << "Ratios are not all integers:" << endl;
			cerr << "Nl = " << Nl << ", Ng = " << p.Ng << ", Nw = " << Nw << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		else {
			Npg = (uint)Nl/p.Ng;
			Npw = (uint)Nl/Nw;
		}
	}

	/*----------------------------------------------------------------------------------------------------------------------------
		4. coordinating files
	----------------------------------------------------------------------------------------------------------------------------*/

	uint loopMin = rank*Npw;
	uint loopMax = (rank+1)*Npw-1;
	
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

	/*----------------------------------------------------------------------------------------------------------------------------
		5. evaluating loop quantitites
	----------------------------------------------------------------------------------------------------------------------------*/

	uint Seed = time(NULL)+rank+2;
	Loop<dim> l(p.K,Seed);
	uint counter = 0;
	uint id;
	number s0, v, w;
	number localSums[2*Nq];

	for (uint j=0; j<Npw; j++) {
		counter++;
		l.load(folder[j]);

		s0 = S0(l);
		v = V0(l);
		w = gsl_sf_cos(p.G*I0(l));
		localSums[0] += s0;
		localSums[1] += s0*s0;
		localSums[2] += v;
		localSums[3] += v*v;
		localSums[4] += w;
		localSums[5] += w*w;
	
		if (counter==Npg) {
			//id = rank*(p.Ng/Nw)+((j+1)/Npg-1);
			MPI_Reduce(&localSums, &sums, 2*Nq, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
			memset(localSums,0,sizeof(localSums));
			counter = 0;
		}
	}

	/*----------------------------------------------------------------------------------------------------------------------------
		6. evaluating errors
	----------------------------------------------------------------------------------------------------------------------------*/

	if (rank==root) {
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
		7. printing results
	----------------------------------------------------------------------------------------------------------------------------*/

		string timenumber = currentDateTime();	
	
		Filename rf = "results/"+timenumber+"loopGroups_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+".dat";
		FILE * ros;
		ros = fopen(((string)rf).c_str(),"w");
		for (uint j=0; j<p.Ng; j++) {
			fprintf(ros,"%12s%5i%5i%8i%8i%10.2g%8i%13.5g%13.5g%13.5g%13.5g%13.5g%13.5g\n",\
					timenumber.c_str(),dim,p.K,Nl,p.Ng,p.G,j,aprxS0_g[j],aprxS0,errorS0,aprxV_g[j],aprxV,errorV);
		}
		fclose(ros);	
		cout << "results printed to " << rf << endl;
	
		rf = "results/loop_dim_"+nts<uint>(dim)+".dat";
		ros = fopen(((string)rf).c_str(),"a");
			fprintf(ros,"%12s%5i%5i%8i%8i%10.2g%13.5g%13.5g%13.5g%13.5g\n",\
					timenumber.c_str(),dim,p.K,Nl,p.Ng,p.G,aprxS0,errorS0,aprxV,errorV);
		fclose(ros);	
		cout << "results printed to " << rf << endl;

		cout << "timenumber: " << timenumber << endl;
		printf("\n");
		printf("%8s%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s\n","dim","Nl","Ng","K","g","S0","varS0",\
			"errorS0","V","varV","errorV");
		printf("%8i%8i%8i%8i%8.2g%12.3g%12.3g%12.3g%12.3g%12.3g%12.3g\n",\
			dim,Nl,p.Ng,p.K,p.G,aprxS0,varianceS0,errorS0,aprxV,varianceV,errorV);
		printf("\n");
	}

}

MPI_Barrier(MPI_COMM_WORLD);
MPI::Finalize();

return 0;
}
