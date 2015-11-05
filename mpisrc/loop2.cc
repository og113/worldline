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
		1 - getting argv
		2 - getting parameters
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

int Nwi, rank, root = 0; // Nw, number of workers
uint Nw;

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &Nwi);
Nw = Nwi;

if (rank==root)
	cout << "starting loop2 with " << Nw << " nodes" << endl;

	
/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

bool verbose = false;
bool circle = false;
string inputsFile = "inputs2";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("circle")==0) circle = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
}

if (rank==root)
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

if (circle)
	circle = (abs(p.G*p.B)>MIN_NUMBER);

MPI_Barrier(MPI_COMM_WORLD);

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
number time_met;

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
			cerr << "Loop2 error: Need Nl = Nw" << endl;
			cerr << "Nl = " << p.Nl << ", Nw = " << Nw << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	/*----------------------------------------------------------------------------------------------------------------------------
		5. coordinating files and data arrays
	----------------------------------------------------------------------------------------------------------------------------*/
	
	// in files
	Filename loadFile = (circle?\
				"data/circle/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_R_"+nts<uint>(p.G*p.B)\
										+"_rank_"+nts<uint>(rank)+".dat":\
				"data/s0+v/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_G_"+nts<uint>(p.G)\
	
									+"_B_"+nts<uint>(p.B)+"_rank_"+nts<uint>(rank)+".dat");
	if (rank==root) {
		cout << "circle " << circle << endl;
		cout << loadFile << endl;
		cout << loadFile.exists() << endl;
	}
	// check if file exists
	if (!loadFile.exists()) {
		loadFile = "data/s0/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_run_"+nts<uint>(rank)+".dat";
		if (!loadFile.exists()) {
			cerr << "Loop2 error: " << loadFile << " doesn't exist" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	if (rank==root) {
		cout << "loading loops from:" << endl;
		cout << loadFile << endl;
	}

	// out files
	string timenumber = currentDateTime();
	Filename loopFile = (circle? (string)loadFile:\
				"data/s0+v/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_G_"+nts<uint>(p.G)\
										+"_B_"+nts<uint>(p.B)+"_rank_"+nts<uint>(rank)+".dat");
	Filename s0File = "data/s0+v/local/"+timenumber+"s0_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+"_G_"+nts<uint>(p.G)\
									+"_B_"+nts<uint>(p.B)+"_rank_"+nts<uint>(rank)+".dat";
	//Filename corrTotalFile = "data/s0+v/vCorrTotal_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+"_G_"+nts<uint>(p.G)\
									+"_B_"+nts<uint>(p.B)+"_rank_"+nts<uint>(rank)+".dat";
	Filename frFile = s0File, vFile = s0File, corrFile = s0File;
	frFile.ID = "fr";
	vFile.ID = "v";
	corrFile.ID = "vCorr";
	
	{
		// local data arrays
		vector<number> s0_data_local(p.Nsw,0.0);
		vector<number> fr_data_local(p.Nsw,0.0);
		vector<number> v_data_local(p.Nsw,0.0);

/*----------------------------------------------------------------------------------------------------------------------------
			6. evaluating loop quantitites
----------------------------------------------------------------------------------------------------------------------------*/

		uint Seed = time(NULL)+rank+2, steps_local = 0, steps;
		Loop<dim> loop(p.K,Seed);
		Metropolis<dim> met(loop,p,++Seed);
		number s0, v, I, fr, lp = 1.0/p.G/p.B;
		
		// timing  metropolis
		clock_t time_run = 0.0;
		if (rank==root)
			time_run = clock();
		
		loop.load(loadFile);
	
		// doing dummy metropolis runs
		met.step(p.Nig*Np);
		met.setSeed(time(NULL)+rank+2);
		
		if (rank==root && verbose)
			printf("%8s%12s%12s%12s%12s\n","sweep","S0","V","I","Fr");
		
		for (uint k=0; k<p.Nsw; k++) {
	
			// metropolis runs per sweep
			if (k==(p.Nsw-1))
				steps_local = met.step(p.Npsw*Np);
			else
				met.step(p.Npsw*Np);
			met.setSeed(time(NULL)+k*1000+rank+2);
		
			s0 = S0(loop);
			I = I0(loop);
			//w = gsl_sf_cos(p.G*I0(loop));
			v = p.G*V1r(loop,p.Epsi);
			v -= (abs(p.Epsi)>MIN_NUMBER? p.G*pi*L(l)/p.Epsi: 0.0);
			fr = (I<lp? 0.0: (-(pi*lp/2.0)*(I-lp/2.0))+pi*I*I/4.0)*gsl_sf_exp(-v);
			//f = (I<lp? -pi*I*I/4.0: -(pi*lp/2.0)*(I-lp/2.0));
		
			s0_data_local[k] = s0;
			fr_data_local[k] = fr;
			v_data_local[k] = v;
			
			if (rank==root && verbose)
				printf("%8i%12.5g%12.5g%12.5g%12.5g\n",k,s0,v,I,fr);
		
		}
		
		// calculating time per successful metropolis step
		MPI_Reduce(&steps_local, &steps, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
		if (rank==root) {
			time_run = clock()-time_run;
			number realtime = time_run/1000000.0;
			time_met = (number)Nw*realtime/steps;
		}
	
/*----------------------------------------------------------------------------------------------------------------------------
			7. printing data
----------------------------------------------------------------------------------------------------------------------------*/
	
		// printing loops
		loop.save(loopFile);
	
		// printing results
		saveVectorBinary(s0File,s0_data_local);
		saveVectorBinary(frFile,fr_data_local);
		saveVectorBinary(vFile,v_data_local);
		if (rank==root)
			cout << "data printed to: " << endl << s0File << endl << frFile << endl << vFile << endl;
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
	vector<number> errorSqrd_local(Nr,0.0), errors;


	// monte carlo data analysis (MCDA)
	MonteCarloData s0MCDA(s0File), frMCDA(frFile), vMCDA(vFile);
	
	// calculating averages
	s0MCDA.calcMeans(avgs_local[0],avgsSqrd_local[0]);
	frMCDA.calcMeans(avgs_local[1],avgsSqrd_local[1]);
	vMCDA.calcMeans(avgs_local[2],avgsSqrd_local[2]);
	
	// calculating correlations
	s0MCDA.calcCorrs(intCorrTime_local[0],expCorrTime_local[0],corrErrorSqrd_local[0]);
	frMCDA.calcCorrs(intCorrTime_local[1],expCorrTime_local[1],corrErrorSqrd_local[1]);
	vMCDA.calcCorrs(intCorrTime_local[2],expCorrTime_local[2],corrErrorSqrd_local[2]);

	// saving correlations
	vMCDA.saveCorrelator(corrFile);
	//vMCDA.saveCorrelatorAppendAscii(corrTotalFile);
	if (rank==root)
			cout << "correlators printed to: " << endl << corrFile << endl;// << corrTotalFile << endl;
	
	// calculating errors
	uint bootstraps = p.Nsw*1e2;
	uint Seed = time(NULL)+rank+2;
	errorSqrd_local[0] = s0MCDA.calcBootstrap(bootstraps,Seed);
	errorSqrd_local[1] = frMCDA.calcBootstrap(bootstraps,Seed);
	errorSqrd_local[2] = vMCDA.calcBootstrap(bootstraps,Seed);
	for (uint k=0; k<Nr; k++) {
		if (abs(errorSqrd_local[k])>MIN_NUMBER)
			weighting_local[k] = 1.0/errorSqrd_local[k];
		else if (abs(avgs_local[k])>MIN_NUMBER)
			cerr << "loop2 error: errorSqrd_local[" << k << "] = 0.0" << endl;
	}
	
	// preparing to combine averages and errors
	for (uint k=0; k<Nr; k++) 
		avgs_local[k] *= weighting_local[k];

	if (rank==root) {
		avgs.resize(Nr,0.0);
		weighting.resize(Nr,0.0);
		errors.resize(Nr,0.0);
		intCorrTime.resize(Nr,0.0);
		expCorrTime.resize(Nr,0.0);
		errors.resize(Nr,0.0);
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
	
		Filename rf = "results/s0+v/loop2_dim_"+nts<uint>(dim)+".dat";
		rf.ID += "Cosmos";
		FILE * ros;
		ros = fopen(((string)rf).c_str(),"a");
		fprintf(ros,"%12s%5i%5i%8i%8i%8.4g%8.4g%8.4g",timenumber.c_str(),dim,p.K,p.Nl,p.Nsw,p.G,p.B,p.Epsi);
		for (uint j=0; j<Nr; j++)
			fprintf(ros,"%13.5g%13.5g%13.5g%13.5g",avgs[j],errors[j],intCorrTime[j],expCorrTime[j]);
		fprintf(ros,"%13.5g",time_met);
		fprintf(ros,"\n");
		fclose(ros);		
		cout << "results printed to:" << endl << rf << endl << endl;
	
		cout << "timenumber: " << timenumber << endl;
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%13s%13s%13s%13s%13s\n","dim","Nl","Nsw","K","G","B",\
				"v","%v","T_int","T_exp","time_met");
		printf("%8i%8i%8i%8i%8.4g%8.4g",dim,p.Nl,p.Nsw,p.K,p.G,p.B);
		printf("%13.5g%13.5g%13.5g%13.5g%13.5g",avgs[2],100.0*errors[2]/abs(avgs[2]),intCorrTime[2],expCorrTime[2],time_met);
		printf("\n\n");
		
	}

}

MPI_Barrier(MPI_COMM_WORLD);
MPI::Finalize();

return 0;
}
