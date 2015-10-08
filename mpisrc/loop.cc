/*----------------------------------------------------------------------------------------------------------------------------
	loop
		program to calculate quantities about worldlines, using gaussian probability distribution
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
		if (id.compare("data")==0 || id.compare("dataChoice")==0) dataChoice = (string)(argv[2*j+2]);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
}

// setting dataChoice to a standard form
if (dataChoice.compare("S0")==0 || dataChoice.compare("s0")==0) {
	dataChoice = "s0";
}
else if (dataChoice.compare("V")==0 || dataChoice.compare("v")==0) {
	dataChoice = "v";
}
else if (dataChoice.compare("Z")==0 || dataChoice.compare("z")==0) {
	dataChoice = "z";
}
else if (dataChoice.compare("W")==0 || dataChoice.compare("w")==0) {
	dataChoice = "w";
}
else if (!dataChoice.empty()) {
	cerr << "dataChoice, " << dataChoice << ", not understood" << endl;
	dataChoice = "";
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. defining basic quantitites	
----------------------------------------------------------------------------------------------------------------------------*/

// parameter loops
uint Npl = 1; // number of parameter loops
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label)) {
	Npl = (pr.Steps)[label-1];
	if (rank==root)
		cout << "looping " << label << " over " << Npl << " steps" << endl;
}

// quantities to calculate sums for
// elements in sums are in pairs of (X,X^2): S0, V, W
// these quantities can be calculated quickly and parallely
uint Nr = 3; // number of different results desired
uint Na = 1; // number of auxialliary quantities to calculate
uint Nct = 1; // number of cross terms to calculate
uint Nq = 2*(Nr+Na)+Nct; // total number of quantities to sum
number *sums = NULL;
number *temp = NULL; // as MPI_SUM doesn't store previous value

// distribution of quantity
number *data = NULL;

/*----------------------------------------------------------------------------------------------------------------------------
	3. starting parameter loop
		- dealing with parameters
		- initializing data arrays
----------------------------------------------------------------------------------------------------------------------------*/

for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pr.toStep(label) && pl>0)
		p.step(pr);
		
	uint Npg = (uint)p.Nl/p.Ng; // Npg, number of loops per group
	uint Npw = (uint)p.Nl/Nw; // Npw, number of loops per worker
	uint Ngpw = (uint)p.Ng/Nw; // Ngpw, number of groups per worker
	
	if (rank==root) {
		// checking relevant ratios are integers
		if (p.Nl%p.Ng!=0 || p.Nl%Nw!=0 || p.Ng%Nw!=0) {
			cerr << "Relevant ratios are not all integers:" << endl;
			cerr << "Nl = " << p.Nl << ", Ng = " << p.Ng << ", Nw = " << Nw << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		// allocating space for data in root
		sums = new number[Nq](); // n.b. () is for initializing to zero
		temp = new number[Nq]();
		if (!dataChoice.empty())
			data = new number[p.Ng]();
	}

	/*----------------------------------------------------------------------------------------------------------------------------
		4. coordinating files and data arrays
	----------------------------------------------------------------------------------------------------------------------------*/
	
	uint loopMin = rank*Npw;
	uint loopMax = (rank+1)*Npw-1;
	
	// constructing folders
	FilenameAttributes faMin, faMax;
	faMin.Directory = "data/loops/s0/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K);
	faMin.Timenumber = "";
	faMax = faMin;
	(faMin.Extras).push_back(StringPair("run",nts<uint>(loopMin)));
	(faMax.Extras).push_back(StringPair("run",nts<uint>(loopMax)));

	Folder folder(faMin,faMax);
	if (folder.size()!=(loopMax-loopMin+1)) {
		cerr << "error for processor " << rank << ":" << endl;
		cerr << "folder.size() = " << folder.size() << ", loopMax-loopMin+1 = " << loopMax-loopMin+1 << endl;
		return 1;
	}
	
	// local data arrays
	number *data_local = NULL;
	if (!dataChoice.empty())
		data_local = new number[Ngpw]();

	/*----------------------------------------------------------------------------------------------------------------------------
		5. evaluating loop quantitites
	----------------------------------------------------------------------------------------------------------------------------*/

	uint Seed = time(NULL)+rank+2;
	Loop<dim> l(p.K,Seed);
	uint counter = 0, id;
	number s0, vz, z, w, gbt = p.G*p.B*p.T;
	number *sums_local = new number[Nq]();

	for (uint j=0; j<Npw; j++) {
		counter++;
		l.load(folder[j]);

		s0 = S0(l);
		vz = p.G;//V0(l);
		z = gsl_sf_exp(-vz);
		vz *= z;
		w = gsl_sf_cos(gbt*I0(l));
		sums_local[0] += s0;
		sums_local[2] += w;
		sums_local[4] += vz;
		sums_local[6] += z;
	
		if (counter==Npg) {
			for (uint k=0; k<(Nr+Na); k++) 
				sums_local[2*k+1] = sums_local[2*k]*sums_local[2*k];
			sums_local[8] = sums_local[4]*sums_local[6];
			
			if (!dataChoice.empty()) {
				id = ((j+1)/Npg-1); // for global id: +rank*(p.Ng/Nw)
				if (dataChoice.compare("s0")==0)
					data_local[id] = sums_local[0]/(number)Npg;
				else if (dataChoice.compare("w")==0) 
					data_local[id] = sums_local[2]/(number)Npg;
				else if (dataChoice.compare("v")==0)
					data_local[id] = sums_local[4]/sums_local[6];
				else if (dataChoice.compare("z")==0)
					data_local[id] = sums_local[4]/(number)Npg;
			}
			
			MPI_Reduce(sums_local, temp, Nq, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
			if (rank==root) {
				for (uint k=0; k<Nq; k++)
					sums[k] += temp[k];
			}
			memset(sums_local,0,Nq*sizeof(number));
			counter = 0;
		}
	}
	
	delete[] sums_local;
	
	// gathering data
	if (!dataChoice.empty())
		MPI_Gather(data_local, Ngpw, MPI_DOUBLE, data, Ngpw, MPI_DOUBLE, root, MPI_COMM_WORLD);

	/*----------------------------------------------------------------------------------------------------------------------------
		6. evaluating errors
	----------------------------------------------------------------------------------------------------------------------------*/

	if (rank==root) {
		vector<number> averages(Nr);
		vector<number> errors(Nr);
		number average2, variance;
		
		// first, results not requiring auxilliary quantities, i.e. s0 and w
		for (uint j=0; j<(Nr-Na); j++) {
			averages[j] = sums[2*j]/(number)p.Nl;
			average2 = sums[2*j+1]/(number)Npg/(number)p.Nl;
			variance = average2-averages[j]*averages[j];
			errors[j] = sqrt(variance/(p.Ng-1.0));
		}
		
		// then the more complicated ones, v from vz and z
		if (Na>0) {
			number v, v2, z, z2, vz;
			number sigma_vv, sigma_vz, sigma_zz;
			
			v = sums[2*(Nr-Na)]/(number)p.Nl;
			v2 = sums[2*(Nr-Na)+1]/(number)Npg/(number)p.Nl;
			z = sums[2*(Nr-Na+1)]/(number)p.Nl;
			z2 = sums[2*(Nr-Na+1)+1]/(number)Npg/(number)p.Nl;
			vz = sums[2*(Nr+Na)]/(number)Npg/(number)p.Nl;
			sigma_vv = v2-v*v;
			sigma_vz = vz-v*z;
			sigma_zz = z2-z*z;
			
			averages[Nr-Na] = v/z;
			errors[Nr-Na] = sqrt((sigma_vv - 2.0*v*sigma_vz/z + v*v*sigma_zz/z/z)/z/z);
			
		}

	/*----------------------------------------------------------------------------------------------------------------------------
		7. printing results
	----------------------------------------------------------------------------------------------------------------------------*/

		string timenumber = currentDateTime();	
	
		Filename rf = "results/loop_dim_"+nts<uint>(dim)+".dat";
		rf.ID += "Cosmos";
		FILE * ros;
		ros = fopen(((string)rf).c_str(),"a");
		fprintf(ros,"%12s%5i%5i%8i%8i%8.2g%8.2g%8.2g",timenumber.c_str(),dim,p.K,p.Nl,p.Ng,p.G,p.B,p.T);
		for (uint j=0; j<Nr; j++)
			fprintf(ros,"%13.5g%13.5g",averages[j],errors[j]);
		fprintf(ros,"\n");
		fclose(ros);
		
		cout << "results printed to " << rf << endl;
		if (!dataChoice.empty()) {
			rf = "data/"+timenumber+"loop_data_"+dataChoice+"_dim_"+nts<uint>(dim)+"_K_"+nts<uint>(p.K)+".dat";
			ros = fopen(((string)rf).c_str(),"w");
			for (uint j=0; j<p.Ng; j++) {
				fprintf(ros,"%12s%5i%5i%8i%8i%8.5g%8.5g%8.5g%8i%13.5g\n",timenumber.c_str(),dim,p.K,p.Nl,p.Ng,p.G,p.B,p.T,j,data[j]);
			}
			fclose(ros);	
			cout << "results printed to " << rf << endl;
		}
	
		cout << "timenumber: " << timenumber << endl;
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s\n","dim","Nl","Ng","K","G","B","T","S0",\
			"%errorS0","W","%errorW","V","%errorV");
		printf("%8i%8i%8i%8i%8.5g%8.5g%8.5g",dim,p.Nl,p.Ng,p.K,p.G,p.B,p.T);
		for (uint j=0; j<Nr; j++)
			printf("%12.4g%12.4g",averages[j],100.0*errors[j]/averages[j]);
		printf("\n\n");
		
		// deleting space for data in root
		delete[] sums;
		sums = NULL;
		delete[] temp;
		temp = NULL;
		delete[] data;
		data = NULL;
	}
	
	delete[] data_local;
	data_local = NULL;
}

MPI_Barrier(MPI_COMM_WORLD);
MPI::Finalize();

return 0;
}
