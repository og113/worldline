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
#include "folder.h"
#include "genloop.h"
#include "parameters.h"
#include "simple.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		0 - getting parameters
		1 - initializing mpi
		2 - getting argv
		2 - defining basic quantitites
		3 - starting parameter loop
		4 - coordinating files
		5 - evaluating loop quantitites
		6 - evaluating errors
		7 - printing results
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

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
if (p.G==0 || p.Nms==0 || p.Nsw==0 || p.Npsw==0 ) {
	cerr << "trivial loop2 run due to parameters: " << endl;
	cerr << p << endl;
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
uint Nq = 2*Nr; // total number of quantities to sum
number *avgs = NULL;
number *weighting = NULL;
number expCorrTime = 0.0, intCorrTime = 0.5;

// full data for quantities
number *data_s0 = NULL, *data_w = NULL, *data_v = NULL;

/*----------------------------------------------------------------------------------------------------------------------------
	3. starting parameter loop
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
			cerr << "Nl = " << Nl ", Nw = " << Nw << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		// allocating space for data in root
		avgs = new number[Nq](); // n.b. () is for initializing to zero
		errorsSqrd = new number[Nq]();
		data_s0 = new number[Nw*p.Nsw]();
		data_w = new number[Nw*p.Nsw]();
		data_v = new number[Nw*p.Nsw]();
	}

	/*----------------------------------------------------------------------------------------------------------------------------
		4. coordinating files and data arrays
	----------------------------------------------------------------------------------------------------------------------------*/
	
	// loop file
	Filename file = "data/gaussian/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_run_"+nts<uint>(rank)+".dat";

	// local data arrays
	number *data_local_s0 = new number[p.Nsw]();
	number *data_local_w = new number[p.Nsw]();
	number *data_local_v = new number[p.Nsw]();
	number *autoCorr = new number[p.Nsw-1](); //note there are one fewer autocorrelations than sweeps

	/*----------------------------------------------------------------------------------------------------------------------------
		5. evaluating loop quantitites
	----------------------------------------------------------------------------------------------------------------------------*/

	uint Seed = time(NULL)+rank+2;
	Loop<dim> loop(p.K,Seed);
	Metropolis<dim> met(loop,p.K,++Seed);
	uint groupCounter = 0, id;
	number s0, vz, z, w;
	number *avgs_local = new number[Nq]();

	loop.load(file);
	
	// doing dummy metropolis runs
	if (abs(p.G)>MIN_NUMBER && p.Nsw>0) {
		met.step(Nig*Np);
		met.setSeed(time(NULL)+rank+2);
	}
	
	for (uint k=0; k<Nsw; k++) {
	
		// metropolis runs per sweep
		if (abs(p.G)>MIN_NUMBER && p.Nms>0) {
			met.step(Npsw*Np);
			met.setSeed(time(NULL)+k*1000+rank+2);
		}
		
		s0 = S0(loop);
		w = gsl_sf_cos(p.G*I0(loop));
		v = V0(loop);
		avgs_local[0] += s0;
		avgs_local[2] += w;
		avgs_local[4] += v;
		for (uint l=0; l<Nr; l++) 
			avgs_local[2*l+1] += avgs_local[2*l]*avgs_local[2*l];
		
		data_local_s0[k] = local[0];
		data_local_w[k] = local[2];
		data_local_v[k] = local[4];
		
	}
	for (uint l=0; l<Nq; l++) 
		avgs_local[l] /= (number)Nsw;
	
	// gathering data - do i really need all this data together?
	MPI_Gather(data_local_s0, Ngpw, MPI_DOUBLE, data_s0, Ngpw, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Gather(data_local_w, Ngpw, MPI_DOUBLE, data_w, Ngpw, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Gather(data_local_v, Ngpw, MPI_DOUBLE, data_v, Ngpw, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
	/*----------------------------------------------------------------------------------------------------------------------------
		6. evaluating errors
	----------------------------------------------------------------------------------------------------------------------------*/

	// calculating (scaled) autocorrelations, using v
	number scaling = avgs_local[5]-avgs_local[4]*avgs_local[4];
	autoCorr[0] = 1.0;
	
	uint expCount = 0;
	bool expBool = true;
	
	// n.b. expensive double sum, could be relegated to another analysis program
	for (uint k=1; k<(Nsw-1); k++) {
		for (uint l=0; l<(Nsw-k); l++)
			autoCorr[k] += data_local_v[l]*data_local_v[l+k];
		autoCorr[k] /= (number)(Nsw-1.0-k)
		autoCorr[k] -= avgs_local[4]*avgs_local[4];
		autoCorr[k] /= scaling;
		
		intCorrTime += autoCorr[k];
		
		if (autoCorr[k]>0 && autoCorr[k]<autoCorr[k-1] && expBool) {
			expCount++;
			expCorrTime += -(number)k/gsl_sf_log(autoCorr[k]);
		}
		else
			expBool = false;
	}
	if (expCount!=0)
		expCorrTime /= (number)expCount;

	// calculating errors locally
	number variance;
	number *weighting_local = new number[Nr]();
	// errors
	for (uint j=0; j<Nr; j++) {
		variance = avgs_local[2*j+1] - avgs_local[2*j]*avgs_local[2*j];
		weighting_local[j] = (number)Nsw/intCorrTime/variance;
		avgs_local[2*j] *= weighting_local[j];
	}
	
	// gathering averages
	MPI_Reduce(avgs_local, avgs, Nq, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	MPI_Reduce(weighting_local, weighting, Nq, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if (rank==root) {
		for (uint k=0; k<Nq; k++)
			avgs[k] /= (number)weighting[k];
	}
	
	delete[] avgs_local;
	avgs_local = NULL;
	delete[] weighting_local;
	weighting_local = NULL;

	/*----------------------------------------------------------------------------------------------------------------------------
		7. printing results
	----------------------------------------------------------------------------------------------------------------------------*/

	if (rank==root) {
		string timenumber = currentDateTime();	
	
		Filename rf = "results/metropolis/loop_dim_"+nts<uint>(dim)+".dat";
		rf.ID += "Office";
		FILE * ros;
		ros = fopen(((string)rf).c_str(),"a");
		fprintf(ros,"%12s%5i%5i%8i%8i%8.2g",timenumber.c_str(),dim,p.K,p.Nl,p.Ng,p.G);
		for (uint j=0; j<Nr; j++)
			fprintf(ros,"%13.5g%13.5g",averages[j],errors[j]);
		fprintf(ros,"\n");
		fclose(ros);
		
		cout << "results printed to " << rf << endl;
		if (!dataChoice.empty()) {
			rf = "data/metropolis/"+timenumber+"loop_data_"+dataChoice+"_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+".dat";
			ros = fopen(((string)rf).c_str(),"w");
			for (uint j=0; j<p.Ng; j++) {
				fprintf(ros,"%12s%5i%5i%8i%8i%8.2g%8i%13.5g\n",timenumber.c_str(),dim,p.K,p.Nl,p.Ng,p.G,j,data[j]);
			}
			fclose(ros);	
			cout << "results printed to " << rf << endl;
		}
	
		cout << "timenumber: " << timenumber << endl;
		printf("\n");
		printf("%8s%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s\n","dim","Nl","Ng","K","G","S0",\
			"%errorS0","W","%errorW","V","%errorV");
		printf("%8i%8i%8i%8i%8.2g",dim,p.Nl,p.Ng,p.K,p.G);
		for (uint j=0; j<Nr; j++)
			printf("%12.4g%12.4g",averages[j],100.0*errors[j]/averages[j]);
		printf("\n\n");
		
		// deleting space for data in root
		delete[] avgs;
		avgs = NULL;
		delete[] weighting;
		weighting = NULL;
		delete[] data_s0;
		data_s0 = NULL;
		delete[] data_w;
		data_w = NULL;
		delete[] data_v;
		data_v = NULL;
	}
	
	delete[] data_local_s0;
	delete[] data_local_w;
	delete[] data_local_v
	delete[] autoCorr;
	data_local_s0 = NULL;
	data_local_w = NULL;
	data_local_v = NULL;
	autoCorr = NULL;
}

MPI_Barrier(MPI_COMM_WORLD);
MPI::Finalize();

return 0;
}
