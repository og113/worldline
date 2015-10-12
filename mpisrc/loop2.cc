/*----------------------------------------------------------------------------------------------------------------------------
	loop2
		program to calculate quantities about worldlines, using metropolis monte-carlo
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
#include <gsl/gsl_sf_log.h>
#include "analysis.h"
#include "folder.h"
#include "genloop.h"
#include "parameters.h"
#include "print.h"
#include "simple.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		0 - initializing mpi
		1 - getting parameters
		2 - getting argv
		3 - defining basic quantitites
		4 - starting parameter loop
		5 - coordinating files
		6 - evaluating loop quantitites
		7 - printing data
		8 - evaluating errors
		9 - printing results
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	0. initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int Nw, rank, root = 0; // Nw, number of workers

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &Nw);

if (rank==root)
	cout << "starting loop with " << Nw << " nodes" << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting parameters
----------------------------------------------------------------------------------------------------------------------------*/

//dimension
#define dim 4

// parameters
ParametersRange pr;
pr.load("inputs");
Parameters p = pr.Min;
if (rank==root) {
	if (p.empty()) {
		cerr << "Parameters empty: nothing in inputs file" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if (abs(p.G)<MIN_NUMBER || p.Nsw==0 || p.Npsw==0 ) {
		cerr << "trivial loop2 run due to parameters: " << endl;
		cerr << p << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
}
MPI_Barrier(MPI_COMM_WORLD);
	
/*----------------------------------------------------------------------------------------------------------------------------
	2. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

// data to print
string dataChoice = "";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("data")==0 || id.compare("dataChoice")==0);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. defining basic quantitites	
----------------------------------------------------------------------------------------------------------------------------*/

// parameter loops
uint Npl = 1; // number of parameter loops
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label) && label!=2) {
	Npl = (pr.Steps)[label-1];
	if (rank==root)
		cout << "looping " << label << " over " << Npl << " steps" << endl;
}

// quantities to calculate sums for
// elements in sums are in pairs of (X,X^2): S0, V, W
// these quantities can be calculated quickly and parallely
uint Nr = 3; // number of different results desired
number mets_local, mets;

/*----------------------------------------------------------------------------------------------------------------------------
	4. starting parameter loop
		- dealing with parameters
		- initializing data arrays
----------------------------------------------------------------------------------------------------------------------------*/

for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pr.toStep(label) && pl>0)
		p.step(pr);
		
	uint Np = pow(2,p.K); // Np, number of points per loop
	
	if (rank==root) {
		// checking Nl==Nw
		if (p.Nl!=Nw) {
			cerr << "Need Nl = Nw for loop2" << endl;
			cerr << "Nl = " << p.Nl << ", Nw = " << Nw << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	/*----------------------------------------------------------------------------------------------------------------------------
		5. coordinating files and data arrays
	----------------------------------------------------------------------------------------------------------------------------*/
	
	// in files
	Filename loadFile = "data/gaussian/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_run_"+nts<uint>(rank)+".dat";
	
	// out files
	Filename loopFile = "results/metropolis/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_B_"+nts<uint>(p.B)\
									+"_G_"+nts<uint>(p.G)+"_rank_"+nts<uint>(rank)+".dat";
	Filename s0File = "results/metropolis/s0_dim_"+nts<uint>(dim)+"K_"+nts<uint>(p.K)+"/s0_B_"+nts<uint>(p.B)\				
									+"_G_"+nts<uint>(p.G)+"_rank_"+nts<uint>(rank)+".dat";
	Filename wFile = s0File, vFile = s0File;
	wFile.ID = "w";
	vFile.ID = "v";
	
	{
		// local data arrays
		vector<number> s0_data_local(p.Nsw,0.0);
		vector<number> w_data_local(p.Nsw,0.0);
		vector<number> v_data_local(p.Nsw,0.0);

/*----------------------------------------------------------------------------------------------------------------------------
			6. evaluating loop quantitites
----------------------------------------------------------------------------------------------------------------------------*/

		uint Seed = time(NULL)+rank+2;
		Loop<dim> loop(p.K,Seed);
		Metropolis<dim> met(loop,p,++Seed);
		number s0, v, w;

		loop.load(loadFile);
	
		// doing dummy metropolis runs
		mets_local = met.step(p.Nig*Np);
		met.setSeed(time(NULL)+rank+2);
	
		for (uint k=0; k<p.Nsw; k++) {
	
			// metropolis runs per sweep
			mets_local += met.step(p.Npsw*Np);
			met.setSeed(time(NULL)+k*1000+rank+2);
		
			s0 = S0(loop);
			w = gsl_sf_cos(p.G*I0(loop));
			v = V0(loop);
		
			s0_data_local[k] = s0;
			w_data_local[k] = w;
			v_data_local[k] = v;
		
		}
		
		// calculating mets, average number of metropolis runs per accepted step
		mets_local /= (1.0+p.Nsw);
		MPI_Reduce(&mets_local, &mets, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
		if (rank==root)
			mets /= (number)Nw;
	
/*----------------------------------------------------------------------------------------------------------------------------
			7. printing data
----------------------------------------------------------------------------------------------------------------------------*/
	
		// printing loops
		loop.save(loopFile);
	
		// printing results
		saveVectorBinary(s0File,s0_data_local);
		saveVectorBinary(wFile,w_data_local);
		saveVectorBinary(vFile,v_data_local);
	
	}
	
	/*----------------------------------------------------------------------------------------------------------------------------
		8. loading results and evaluating errors
	----------------------------------------------------------------------------------------------------------------------------*/
	
	// quantities to calculate
	vector<number> avgs, avgs_local(Nr,0.0), avgsSqrd_local(Nr,0.0);
	vector<number> weighting, weighting_local(Nr,0.0);
	vector<number> intCorrTime_local(Nr,0.0), intCorrTime;
	vector<number> expCorrTime_local(Nr,0.0), expCorrTime;
	vector<number> corrErrorSqrd_local(Nr,0.0), corrErrorSqrd;
	vector<number> errors;

	// monte carlo data analysis (MCDA)
	MonteCarloData s0MCDA(s0File), wMCDA(wFile), vMCDA(vFile);
	
	// calculating averages
	s0MCDA.calcMeans(avgs_local[0],avgsSqrd_local[0]);
	wMCDA.calcMeans(avgs_local[1],avgsSqrd_local[1]);
	vMCDA.calcMeans(avgs_local[2],avgsSqrd_local[2]);
	
	// calculating correlations
	s0MCDA.calcCorrs(intCorrTime[0],expCorrTime[0],corrErrorSqrd[0]);
	wMCDA.calcCorrs(intCorrTime[1],expCorrTime[1],corrErrorSqrd[1]);
	vMCDA.calcCorrs(intCorrTime[2],expCorrTime[2],corrErrorSqrd[2]);
	
	// calculating errors
	for (uint k=0; k<Nr; k++) 
		weighting_local[k] = 1.0/corrErrorSqrd[k]; // n.b. this will change to jacknife or bootstrap once written
	
	// preparing to combine averages and errors
	for (uint k=0; k<Nr; k++) 
		avgs_local[k] *= weighting_local[k];

	if (rank==root) {
		avgs.resize(Nr,0.0);
		weighting.resize(Nr,0.0);
		errors.resize(Nr,0.0);
		intCorrTime.resize(Nr,0.0);
		expCorrTime.resize(Nr,0.0);
	}
	
	// gathering averages over workers
	MPI_Reduce(&avgs_local[0], &avgs[0], Nr, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	MPI_Reduce(&weighting_local[0], &weighting[0], Nr, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	MPI_Reduce(&intCorrTime_local[0], &intCorrTime[0], Nr, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	MPI_Reduce(&expCorrTime_local[0], &expCorrTime[0], Nr, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if (rank==root) {
		for (uint k=0; k<Nr; k++) {
			avgs[k] /= (number)weighting[k];
			errors[k] = (Nw==1? sqrt(1.0/weighting[k]) : sqrt((number)(Nw-1.0)/weighting[k]) );
			intCorrTime[k] /= (number)Nw;
			expCorrTime[k] /= (number)Nw;
		}
	}

	/*----------------------------------------------------------------------------------------------------------------------------
		9. printing results
	----------------------------------------------------------------------------------------------------------------------------*/

	if (rank==root) {
		string timenumber = currentDateTime();	
	
		Filename rf = "results/metropolis/loop2_dim_"+nts<uint>(dim)+".dat";
		rf.ID += "Office";
		FILE * ros;
		ros = fopen(((string)rf).c_str(),"a");
		fprintf(ros,"%12s%5i%5i%8i%8i%8.4g%8.4g",timenumber.c_str(),dim,p.K,p.Nl,p.Ng,p.G,p.B);
		for (uint j=0; j<Nr; j++)
			fprintf(ros,"%13.5g%13.5g%13.5g%13.5g",avgs[j],errors[j],intCorrTime[j],expCorrTime[j]);
		fprintf(ros,"\n");
		fclose(ros);		
		cout << "results printed to " << rf << endl;
	
		cout << "timenumber: " << timenumber << endl;
		printf("\n");
		printf("%8s%8s%8s%8s%8s%13s%13s%13s%13s\n","dim","Nl","Ng","K","G",\
				"W","%errorW","T_int","T_exp");
		printf("%8i%8i%8i%8i%8.4g%8.4g",dim,p.Nl,p.Ng,p.K,p.G,p.B);
		printf("%13.5g%13.5g%13.5g%13.5g",avgs[1],100.0*errors[1]/avgs[1],intCorrTime[1],expCorrTime[1]);
		printf("\n\n");
		
	}
	
}

MPI_Barrier(MPI_COMM_WORLD);
MPI::Finalize();

return 0;
}
