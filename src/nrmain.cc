/*
	main program to run serial n-r method to solve classical monopole worldline equations
*/

#include "check.h"
#include "folder.h"
#include "genloop.h"
#include "nrloop.h"
#include "print.h"
#include "simple.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - print options enum
		1 - argv, parameters etc
		2 - beginning parameter loop
		3 - defining quantities
		4 - beginning nr loop
		5 - assigning dS, ddS etc
		6 - some checks
		7 - print early 1
		8 - solve for delta
		9 - print early 2
		10 - convergence
		11 - print output
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------------------
	0 - print options enum
----------------------------------------------------------------------------------------------------------------------------*/

struct PrintOptions {
	enum Option { none, all, x, mds};
};

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1 - argv, parameters etc
----------------------------------------------------------------------------------------------------------------------------*/

// inputs file
bool verbose = true;
bool circle = true;
string printOpts = "";
string inputsFile = "inputs3";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("circle")==0) circle = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("print")==0) printOpts = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

cout << "using inputs file " << inputsFile << endl;

PrintOptions::Option po = PrintOptions::none;
if (!printOpts.empty()) {
	if (printOpts.compare("all")==0)
		po = PrintOptions::all;
	else if (printOpts.compare("x")==0)
		po = PrintOptions::x;
	else if (printOpts.compare("mds")==0)
		po = PrintOptions::mds;
	else
		cerr << "print options not understood: " << printOpts << endl;
}
//dimension
#define dim 4

// parameters
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}

// timenumber
string timenumber = currentDateTime();

/*----------------------------------------------------------------------------------------------------------------------------
	2 - beginning parameter loop
----------------------------------------------------------------------------------------------------------------------------*/

// parameter loops
uint Npl = 1; // number of parameter loops
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label) && label==Parameters::k) {
	Npl = (pr.Steps)[label-1];
	cout << "looping " << label << " over " << Npl << " steps" << endl;
}

// starting loop
for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pr.toStep(label) && pl>0)
		p.step(pr);
	else if (pr.toStep(label) && label==Parameters::nl)
		p.Nl = (pr.Max).Nl;
	uint N = pow(2,p.K);
	uint zm = dim;
	uint NT = N*dim+zm;
		
/*----------------------------------------------------------------------------------------------------------------------------
	3 - defining quantitites
----------------------------------------------------------------------------------------------------------------------------*/

	//defining a time and starting the clock
	clock_t time;
	time = clock();
	
	// defining checks
	Check checkSol("solution",1.0e-7);
	Check checkSolMax("solution max",1.0e-6);
	Check checkDelta("delta",1.0);
	Check checkSm("smoothness",1.0);
	Check checkInv("inverse",1.0e-16*NT*NT);
	
	// defining scalar quantities
	number len, i0, s, sm;
	
	// defining vector and matrix quantities
	vec x(N*dim);
	vec mds(NT);
	mat dds(NT,NT);
	
	// defining xLoop
	uint Seed = time;
	Loop<dim> xLoop(p.K,Seed);
	
	// x file
	Filename loadFile = (circle?\
				"data/circle/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_R_"+nts<uint>(p.G*p.B)\
										+"_rank_0.dat":\
				"data/s0/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(p.K)+"/loop_run_0.dat");
	// check if file exists
	if (!loadFile.exists()) {
		cerr << "Loop2 error: " << loadFile << " doesn't exist" << endl;
		return 1;
	}
	cout << "loading loops from:" << endl;
	cout << loadFile << endl;
	
	// loading x
	loadVectorBinary(loadFile,x);
	x.conservativeResize(NT);
	for (uint mu=0; mu<zm; mu++)
		x[N*dim+mu] = 0.1;
	
	//defining some quantities used to stop n-r loop
	uint runsCount = 0;
	uint minRuns = 1;
	uint maxRuns = 100;
	
/*----------------------------------------------------------------------------------------------------------------------------
	4 - beginning newton-raphson loop
----------------------------------------------------------------------------------------------------------------------------*/

	while ((!checkSol.good() || !checkSolMax.good() || runsCount<minRuns) && runsCount<maxRuns)  {
		runsCount++;
		
		// initializing to zero
		mds = Eigen::VectorXd::Zero(NT);
		dds = Eigen::MatrixXd::Zero(NT,NT);
		
		// loading x to xLoop - messier than it should be (should work with either a vec or a Loop really)
		vectorToLoop(x,xLoop);	

/*----------------------------------------------------------------------------------------------------------------------------
	5 - assigning mdx and ddx
----------------------------------------------------------------------------------------------------------------------------*/
		
		// some simple scalars
		uint j, k, mu, nu;
		number gb = p.G*p.B;
		
		// bulk
		for (j=0; j<N; j++) {
			for (mu=0; mu<dim; mu++) {
			
				mdS0_nr(j,mu,xLoop,1.0,mds);
				//mdL_nr(j,mu,xLoop,1.0,mds);
				mdI_nr(j,mu,xLoop,-gb,mds);
				
				for (k=0; k<N; k++) {
					for (nu=0; nu<dim; nu++) { // doing a full second loop for v
					
						ddS0_nr(j,mu,k,nu,xLoop,1.0,dds);
						//ddL_nr(j,mu,k,nu,xLoop,1.0,dds); // check symmetry in (j,mu)<->(k,nu)
						ddI_nr(j,mu,k,nu,xLoop,-gb,dds);
						
					}
				}
			}		
		}
		
		// lagrange multiplier terms
		for (j=0; j<N; j++) {
			for (mu=0; mu<dim; mu++) {
			
				mds(N*dim+mu) -= x[j*dim+mu];
				mds(j*dim+mu) -= x[N*dim+mu];
				
				dds(j*dim+mu,N*dim+mu) += 1.0;
				dds(N*dim+mu,j*dim+mu) += 1.0;
				
			}
		}
		
		// assigning scalar quantities - messier than it should be, better to be in sync with the previous loops
		len = L(xLoop);
		i0 = -gb*I0(xLoop);
		s = len + i0;
	
/*----------------------------------------------------------------------------------------------------------------------------
	6 - some checks
----------------------------------------------------------------------------------------------------------------------------*/	

		// smoothness
		sm = Sm(xLoop);
		checkSm.add(sm);
	
/*----------------------------------------------------------------------------------------------------------------------------
	7 - print early 1
----------------------------------------------------------------------------------------------------------------------------*/	
		
		if (po!=PrintOptions::none) {
			Filename early = "data/temp/"+timenumber+"xEarly1_K_"+nts(p.K)+"_G_"+nts(p.G)+"_B_"+nts(p.B)+"_run_"+nts(runsCount)+".dat";
			if (po==PrintOptions::x || po==PrintOptions::all) {
				printAsLoop(early,dim,x);
				printf("%12s%50s\n","x:",((string)early).c_str());
			}
			if (po==PrintOptions::mds || po==PrintOptions::all) {
				early.ID = "mdsEarly1";
				printAsLoop(early,dim,mds);
				printf("%12s%50s\n","mds:",((string)early).c_str());
			}
			
		}
		
/*----------------------------------------------------------------------------------------------------------------------------
	8 - solve for delta
----------------------------------------------------------------------------------------------------------------------------*/	
				
		// initializing delta
		vec delta(NT);
		
		// solving for delta = DDS^{-1}*mdS
		delta = dds.partialPivLu().solve(mds);
		
		
		//independent check on whether calculation worked		
		number invError = (dds*delta - mds).norm();
		checkInv.add(invError);
		checkInv.checkMessage();
		if (!checkInv.good())
			return 1;

		//assigning values to x
		x += delta;
	
/*----------------------------------------------------------------------------------------------------------------------------
	9 - printing early 2
----------------------------------------------------------------------------------------------------------------------------*/	

		if (po!=PrintOptions::none) {
			Filename early = "data/temp/"+timenumber+"xEarly2_K_"+nts(p.K)+"_G_"+nts(p.G)+"_B_"+nts(p.B)+"_run_"+nts(runsCount)+".dat";
			if (po==PrintOptions::x || po==PrintOptions::all) {
				printAsLoop(early,dim,x);
				printf("%12s%50s\n","x:",((string)early).c_str());
			}
		}

/*----------------------------------------------------------------------------------------------------------------------------
	10 - convergence checks
----------------------------------------------------------------------------------------------------------------------------*/	
	
		// calculating norms etc
		number normx = x.norm();
		number normmds = mds.norm();
		number normdelta = delta.norm();
		number maxmds = mds.maxCoeff();
		number minmds = mds.minCoeff();
		if (-minmds>maxmds) maxmds = -minmds;
		number maxx = x.maxCoeff();
		number minx = x.minCoeff();
		if (-minx>maxx) maxx = -minx;

		// adding to checks
		checkSol.add(normmds/normx);
		checkSolMax.add(maxmds/maxx);
		checkDelta.add(normdelta/normx);
		
		// checking delta
		checkDelta.checkMessage();
		if (!checkDelta.good()) {
			break;
		}
	
		//printing tests to see convergence
		if (verbose) {
			if (runsCount==1) {
				printf("%5s%5s%12s%12s%12s%12s%12s%12s%12s\n","pl","run","len","i0","s","sol","solM","delta","Sm");
			}
			printf("%5i%5i%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g\n",pl,runsCount,len,i0,s,checkSol.back(),\
				checkSolMax.back(),checkDelta.back(),checkSm.back());
		}
	
	}
	
/*----------------------------------------------------------------------------------------------------------------------------
	11 - printing results
----------------------------------------------------------------------------------------------------------------------------*/	

	//stopping clock
	time = clock() - time;
	number realtime = time/1000000.0;
	
	// printing results to terminal
	printf("\n");
	printf("%8s%8s%8s%8s%8s%8s%14s%14s%14s\n","runs","time","K","G","B","a","len",\
		"i0","s");
	printf("%8i%8g%8i%8g%8g%8g%14.5g%14.5g%14.5g\n",\
		runsCount,realtime,p.K,p.G,p.B,p.Epsi,len,i0,s);
	printf("\n");
	
	
	if (checkDelta.good()) {
	
		// printing results to file	
		string resFile = "results/nr/nrmain.dat";
		FILE* ros;
		ros = fopen(resFile.c_str(),"a");
		fprintf(ros,"%12s%8i%8i%8g%8g%8g%16.6g%16.6g%16.6g%12.4g\n",\
					timenumber.c_str(),pl,p.K,p.G,p.B,p.Epsi,len,i0,s,checkSol.back());
		fclose(ros);
		printf("%12s%50s\n","results:",resFile.c_str());
		
		// printing loop to file
		string loopRes = "data/nr/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_G_"+nts(p.G)+"_B_"+nts(p.B)+".dat";
		saveVectorBinary(loopRes,x);
		printf("%12s%50s\n","x:",loopRes.c_str());
		
	}	

	// printing extras to ascii files
	if (po!=PrintOptions::none) {
		Filename file = "data/temp/"+timenumber+"x_K_"+nts(p.K)+"_G_"+nts(p.G)+"_B_"+nts(p.B)+"_run_"+nts(runsCount)+".dat";
		if (po==PrintOptions::x || po==PrintOptions::all) {
			saveVectorAscii(file,x);
			printf("%12s%50s\n","x:",((string)file).c_str());
		}
		else if (po==PrintOptions::mds || po==PrintOptions::all) {
			file.ID = "mds";
			saveVectorAscii(file,mds);
			printf("%12s%50s\n","mds:",((string)file).c_str());
		}
	}

}

return 0;
}
