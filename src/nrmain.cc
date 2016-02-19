/*
	main program to run serial n-r method to solve classical monopole worldline equations
*/

#include <Eigen/Dense>
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
	enum Option { none, all, x, mds, dds, delta, curvature};
};

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1 - argv, parameters etc
----------------------------------------------------------------------------------------------------------------------------*/

// argv options
bool verbose = true;
bool lemon = false;
bool step = true;
bool weak = false;
bool eigen = false;
bool curvature = false;
bool old = true;
bool alltests = false; // doing alltests
string printOpts = "";
string inputsFile = "inputs4";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("lemon")==0) lemon = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("step")==0) step = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("weak")==0) weak = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("eigen")==0) eigen = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("curvature")==0) curvature = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("old")==0) old = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("print")==0) printOpts = (string)argv[2*j+2];
		else if (id.compare("alltests")==0) alltests = (stn<uint>(argv[2*j+2])!=0);
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
	else if (printOpts.compare("dds")==0)
		po = PrintOptions::dds;
	else if (printOpts.compare("delta")==0)
		po = PrintOptions::delta;
	else if (printOpts.compare("curvature")==0) {
		po = PrintOptions::curvature;
		curvature = true;
	}
	else
		cerr << "print options not understood: " << printOpts << endl;
}
//dimension
#define dim 4

// parameters
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min, pold = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file" << endl;
	return 1;
}

// timenumber
string timenumber = currentDateTime();
cout << "timenumber: " << timenumber << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	2 - beginning parameter loop
----------------------------------------------------------------------------------------------------------------------------*/

// parameter loops
uint Npl = pr.totalSteps();
cout << "looping over " << Npl << " steps" << endl;

// starting loop
for (uint pl=0; pl<Npl; pl++) {

	// doing things before parameter step
	Filename loadFile, stepFile;
	
	// stepping parameters
	if (pl>0) {
		p = pr.position(pl);
		pold = pr.neigh(pl);
		stepFile = filenameLoopNR<dim>(pold);
		if (pold.Ng>0)
			(stepFile.Extras).push_back(StringPair("Ng",nts(pold.Ng)));
	}
		
	
	// defining some derived parameters	
	uint N = pow(2,p.K);
	uint zm = dim;
	uint NT = N*dim+zm;
	number M = p.P2, R = abs(1.0/p.G/p.B);
	Point<dim> P;
	P[0] = p.P1;
	P[1] = p.P2;
	if (dim>2) {
		P[2] = p.P3;
		M = p.P3;
		if (dim>3) {
			P[3] = p.P4;
			M = p.P4;
		}
	}

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
	Check checkDX("dx<<a",2.0e-1);
	Check checkICMax("ic_max<<1",5.0e-1);
	Check checkICAvg("ic_avg<<1",2.0e-1);
	Check checkCCMax("cc_max<<1",5.0e-1);
	Check checkKgAMax("a*kg_max<<1",5.0e-1);
	Check checkKgAAvg("a*kg_avg<<1",2.0e-1);
	Check checkKgDxMax("dx*kg_max<<1",5.0e-1);
	Check checkKgDxAvg("dx*kg_avg<<1",2.0e-1);
	Check checkSym("symmetric",1.0e-16*NT*NT);
	Check checkInv("inverse",1.0e-16*NT*NT*NT);
	Check checkNegEigenvalue("negative eigenvalue",0.1);
	Check checkNegEigenvector("negative eigenvector",0.1);
	Check checkGamma("gamma",0.1);
	
	// defining scalar quantities
	number len, i0, s, sm, v, vr, fgamma, gamma0, gamma1, ic_max, ic_avg, cc_max, kg_max, kg_avg;
	
	// defining vector and matrix quantities
	vec x(N*dim);
	vec mds(NT);
	mat dds(NT,NT);
	
	// curvature
	vec sc_vec;
	vec kg_vec;
	
	// defining xLoop
	Loop<dim> xLoop(p.K,0);
	
	// x file
	if (pl==0 || !step) {
		if (lemon)
			loadFile = "data/lemon/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R)+"_M_"+nts(M)+"_rank_0.dat";
		else {
			loadFile = filenameLoopNR<dim>(p);
			if (p.Ng>0)
				(loadFile.Extras).push_back(StringPair("Ng",nts(p.Ng)));
			if (!loadFile.exists())
				loadFile = filenameLoopNR<dim>(p);
			if (!loadFile.exists()) {
				if (pl>0)
					loadFile = stepFile;
					if (!loadFile.exists())
						loadFile = filenameLoopNR<dim>(pold);
				if (!loadFile.exists() && old)
					loadFile = "data/nr/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_G_"+nts(p.G)+"_B_"+nts(p.B)+"_M_"+nts(p.P4)\
		+"_a_"+nts(p.Epsi)+".dat";
				if (!loadFile.exists())
					loadFile = "data/circle/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R)+"_rank_0.dat";
			}
		}
			
	}
	else {
		loadFile = stepFile;
		lemon = false;
	}
	// check if file exists
	if (!loadFile.exists()) {
		cerr << "nrmain error: " << loadFile << " doesn't exist" << endl;
		return 1;
	}
	cout << "loading loops from:" << endl;
	cout << loadFile << endl;
	
	// loading x
	loadVectorBinary(loadFile,x);
	if (pold.K!=p.K && loadFile==stepFile) {
		Loop<dim> interpOld(pold.K,0), interpNew(p.K,0);
		vectorToLoop(x,interpOld);
		interpolate(interpOld,interpNew);
		loopToVector(interpNew,x);
	}
	if (x.size()!=NT) {
		x.conservativeResize(NT);
		for (uint mu=0; mu<zm; mu++)
			x[N*dim+mu] = 1.0e-3;
	}
	
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
		len = 0.0, i0 = 0.0, v = 0.0, fgamma = 0.0, gamma0 = 0.0, gamma1 = 0.0;//, s0 = 0.0;
		ic_max = 0.0, ic_avg = 0.0, cc_max = 0.0, kg_max = 0.0, kg_avg = 0.0;
		if (curvature) {
			sc_vec = Eigen::VectorXd::Zero(N);
			kg_vec = Eigen::VectorXd::Zero(N);
		}
				
		// loading x to xLoop - messier than it should be (should work with either a vec or a Loop really)
		vectorToLoop(x,xLoop);

/*----------------------------------------------------------------------------------------------------------------------------
	5 - assigning mdx and ddx
----------------------------------------------------------------------------------------------------------------------------*/
		
		// some simple scalars
		uint j, k, mu, nu;
		number gb = p.G*p.B;
		number g = p.G*p.G/8.0/PI/PI;
		number sqrt4s0 = 2.0*sqrt(S0(xLoop));
		number cusp_scale = g*log(p.Mu/p.Epsi);
		number ic_scale = 1.0;
		number cc_scale = 1.0;
		number kg_scale = 1.0;
		
		Point<dim> P0;
		
		// bulk
		for (j=0; j<N; j++) {
		
			//S0(j, xLoop, s0norm, s0norm);
			L		(j, xLoop, 1.0, len);
			I0		(j, xLoop, -gb, i0);
			
			//curvatures
			if (curvature) {
					InlineCurvatureMax(j, xLoop, ic_scale, sc_vec[j]);
					KGMaxPlane(j, xLoop, kg_scale, kg_vec[j]);
			}
			if (!(P^=P0)) {
				InlineCurvatureMax(j, xLoop, 0, N/2-1, ic_scale, ic_max);
				InlineCurvatureAvg(j, xLoop, 0, N/2-1, ic_scale, ic_avg);
				KGMaxPlane(j, xLoop, 0, N/2-1, kg_scale, kg_max);
				KGAvgPlane(j, xLoop, 0, N/2-1, kg_scale, kg_avg);
			}
			else {
				InlineCurvatureMax(j, xLoop, ic_scale, ic_max);
				InlineCurvatureAvg(j, xLoop, ic_scale, ic_avg);
				CuspCurvatureMax(j,xLoop,cc_scale,ic_max);
				KGMaxPlane(j, xLoop, kg_scale, kg_max);
				KGAvgPlane(j, xLoop, kg_scale, kg_avg);
			}
		
			for (mu=0; mu<dim; mu++) {
			
				// free particle
				mdsqrtS0_nr(j,mu,xLoop,sqrt4s0,1.0,mds);
				
				// external field
				mdI_nr(j,mu,xLoop,-gb,mds);
				
				// dynamical field self-energy regularisation
				if (!weak) mdL_nr(j,mu,xLoop,-g*PI/p.Epsi,mds);
				
				for (k=0; k<N; k++) {
				
					if (mu==0)
						V1r(j, k, xLoop, p.Epsi, g, v);
					
					// dynamical field
					if (!weak) mdV1r_nr(j, mu, k, xLoop, p.Epsi, g, mds);
				
					for (nu=0; nu<dim; nu++) {
					
						// free particle
						ddsqrtS0_nr(j,mu,k,nu,xLoop,sqrt4s0,1.0,dds);
						
						// external field
						ddI_nr(j,mu,k,nu,xLoop,-gb,dds);
						
						// dynamical field	
						if (!weak) ddV1r_nr(j, mu, k, nu, xLoop, p.Epsi, g, dds);
						
						// dynamical field self-energy regularisation
						if (!weak) ddL_nr(j,mu,k,nu,xLoop,-g*PI/p.Epsi,dds);
						
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
		
		if (!(P^=P0)) {
			// external momenta
			mdPX_nr(xLoop,N/2-1,P,-1.0,mds);
			mdPX_nr(xLoop,0,P,1.0,mds);
		
			// dynamical field cusp regularisation
			uint range = 2*p.Ng+1;
			uint left = N/2-1 - (int)p.Ng;
			uint right = 0 + (int)(p.Ng>0)*(N-p.Ng);
			
			for (uint k=0; k<range; k++) {
				Gamma(left, xLoop, 1.0, gamma1);
				Gamma(right, xLoop, 1.0, gamma0);
				FGamma(left, xLoop, -cusp_scale, fgamma);
				FGamma(right, xLoop, -cusp_scale, fgamma);
				CuspCurvatureMax(left,xLoop,cc_scale,cc_max);
				CuspCurvatureMax(right,xLoop,cc_scale,cc_max);
				if (!weak) {
					mdFGamma_nr(xLoop,left,-cusp_scale,mds);
					mdFGamma_nr(xLoop,right,-cusp_scale,mds);
					ddFGamma_nr(xLoop,left,-cusp_scale,dds);
					ddFGamma_nr(xLoop,right,-cusp_scale,dds);
				}
				left = (left==(N-1)? 0: left+1);
				right = (right==(N-1)? 0: right+1);
			}
		}
		
		// assigning scalar quantities
		vr = v;
		vr -= (abs(p.Epsi)>MIN_NUMBER? g*PI*len/p.Epsi : 0.0);
		vr -= (!(P^=P0)? cusp_scale*fgamma : 0.0);
		s = sqrt4s0 + i0;
		if (!weak) s += vr;
		if (!(P^=P0))
			s -= Dot(xLoop[N/2-1]-xLoop[0],P);
	
/*----------------------------------------------------------------------------------------------------------------------------
	6 - some checks
----------------------------------------------------------------------------------------------------------------------------*/	

		// smoothness
		sm = Sm(xLoop,0,N/2-1);
		checkSm.add(sm);
		
		// check dx<<a, using the average
		number dx = len/(number)N;
		checkDX.add(dx/p.Epsi);
		
		// check ic<<1, a*kg<<1, dx*kg<<1, cc<<1
		checkICMax.add(ic_max);
		checkICAvg.add(ic_avg);
		checkCCMax.add(cc_max);
		checkKgAMax.add(p.Epsi*kg_max);
		checkKgAAvg.add(p.Epsi*kg_avg);
		checkKgDxMax.add(dx*kg_max);
		checkKgDxAvg.add(dx*kg_avg);
				
		if (alltests) {
			// checking if dds is symmetric
			mat dds_asym(NT,NT);
			dds_asym = dds.transpose();
			dds_asym -= dds;
			dds_asym *= 0.5;
			number asym = dds_asym.norm();
			checkSym.add(asym);
			checkSym.checkMessage();
			
			// checking value of negative eigenvalue, assuming x is an eigenvector
			number analyticNegEigenvalue = -2.0*PI*p.G*p.B;
			number xnorm = x.squaredNorm();
			number negEigenvalue = ((dds*x).dot(x)/xnorm)*(number)xLoop.size();
			number negEigenvalueTest = abs(negEigenvalue-analyticNegEigenvalue)/abs(analyticNegEigenvalue);
			checkNegEigenvalue.add(negEigenvalueTest);
			checkNegEigenvalue.checkMessage();
			
			// checking x is an eigenvector (case for circle in continuum limit)
			negEigenvalue /= (number)xLoop.size();
			number negEigenvectorTest = (dds*x-negEigenvalue*x).squaredNorm()/xnorm;
			checkNegEigenvector.add(negEigenvectorTest);
			checkNegEigenvector.checkMessage();
			
			if (!checkNegEigenvalue.good() || !checkNegEigenvector.good()) {
				cerr << "negative eigenvalue = " << negEigenvalue << endl;
				cerr << "analytic result     = " << analyticNegEigenvalue << endl;
				cerr << "eigenvector test    = " << negEigenvectorTest << endl;
			}
			
			// checking if angle gamma agrees with weak coupling result
			number gamma_ratio = 0.0;
			if (abs(M)>MIN_NUMBER && abs(M)<2.0) {
				number gamma_weak = PI-2.0*asin(sqrt(1.0-pow(M,2)/4.0));
				gamma_ratio = (gamma0+gamma1)/2.0/gamma_weak;
			}
			checkGamma.add(gamma_ratio);
			checkGamma.checkMessage();
		}		
	
/*----------------------------------------------------------------------------------------------------------------------------
	7 - print early 1
----------------------------------------------------------------------------------------------------------------------------*/	
		
		if (po!=PrintOptions::none) {
			Filename early = "data/temp/"+timenumber+"xEarly1_K_"+nts(p.K)+"_G_"+nts(p.G)+"_B_"+nts(p.B)+"_M_"+nts(M)\
						+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+"_run_"+nts(runsCount)+".dat";
			if (po==PrintOptions::x || po==PrintOptions::all) {
				printAsLoop(early,dim,x,N*dim);
				printf("%12s%50s\n","x:",((string)early).c_str());
			}
			if (po==PrintOptions::mds || po==PrintOptions::all) {
				early.ID = "mdsEarly1";
				printAsLoop(early,dim,mds,N*dim);
				printf("%12s%50s\n","mds:",((string)early).c_str());
			}
			if (po==PrintOptions::dds || po==PrintOptions::all) {
				early.ID = "ddsEarly1";
				saveMatrixAscii(early,dds);
				printf("%12s%50s\n","dds:",((string)early).c_str());
			}
			if (po==PrintOptions::curvature || po==PrintOptions::all) {
				early.ID = "curvatureEarly1";
				printAsLoop(early,dim,x,N*dim);
				saveVectorAsciiAppend(early,sc_vec);
				saveVectorAsciiAppend(early,kg_vec);
				printf("%12s%50s\n","curvature:",((string)early).c_str());
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
		if (!checkInv.good()) {
			number x_end = 0.0;
			for (j=0; j<dim; j++)
				x_end += pow(x[N*dim+j],2);
			x_end = sqrt(x_end);
			cerr << endl << "x.norm()         " << x.norm() << endl;
			cerr << "x_end.norm()     " << x_end << endl;
			cerr << "mds.norm()       " << mds.norm() << endl;
			cerr << "mds.maxCoeff()   " << mds.maxCoeff() << endl;
			cerr << "mds.minCoeff()   " << mds.minCoeff() << endl;
			cerr << "delta.norm()     " << delta.norm() << endl;
			cerr << endl << "dds info:" << endl;
			cerr << "dds.determinant():" << dds.determinant() << endl;
			cerr << "dds.sum():        " << dds.sum()       << endl;
			cerr << "dds.prod():       " << dds.prod()      << endl;
			cerr << "dds.mean():       " << dds.mean()      << endl;
			cerr << "dds.minCoeff():   " << dds.minCoeff()  << endl;
			cerr << "dds.maxCoeff():   " << dds.maxCoeff()  << endl;
			cerr << "dds.trace():      " << dds.trace()     << endl;
			cerr << "dds.norm():       " << dds.norm() << endl;
			return 1;
		}

		//assigning values to x
		x += delta;
	
/*----------------------------------------------------------------------------------------------------------------------------
	9 - printing early 2
----------------------------------------------------------------------------------------------------------------------------*/	

		if (po!=PrintOptions::none) {
			Filename early = "data/temp/"+timenumber+"deltaEarly2_K_"+nts(p.K)+"_G_"+nts(p.G)+"_B_"+nts(p.B)+"_M_"+nts(M)\
							+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+"_run_"+nts(runsCount)+".dat";
			if (po==PrintOptions::delta || po==PrintOptions::all) {
				printAsLoop(early,dim,delta,N*dim);
				printf("%12s%50s\n","delta:",((string)early).c_str());
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
			number x_end = 0.0;
			for (j=0; j<dim; j++)
				x_end += pow(x[N*dim+j],2);
			x_end = sqrt(x_end);
			cerr << endl << "x.norm()         " << x.norm() << endl;
			cerr << "x_end.norm()     " << x_end << endl;
			cerr << "mds.norm()       " << mds.norm() << endl;
			cerr << "mds.maxCoeff()   " << mds.maxCoeff() << endl;
			cerr << "mds.minCoeff()   " << mds.minCoeff() << endl;
			cerr << "delta.norm()     " << delta.norm() << endl;
			cerr << endl << "dds info:" << endl;
			cerr << "dds.determinant():" << dds.determinant() << endl;
			cerr << "dds.sum():        " << dds.sum()       << endl;
			cerr << "dds.prod():       " << dds.prod()      << endl;
			cerr << "dds.mean():       " << dds.mean()      << endl;
			cerr << "dds.minCoeff():   " << dds.minCoeff()  << endl;
			cerr << "dds.maxCoeff():   " << dds.maxCoeff()  << endl;
			cerr << "dds.trace():      " << dds.trace()     << endl;
			cerr << "dds.norm():       " << dds.norm() << endl;
			break;
		}
	
		//printing tests to see convergence
		if (verbose) {
			if (runsCount==1) {
				printf("%5s%5s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n","pl","run","len","i0","s","sol","solM","delta","Sm","dx*kg_max","ic_max","cc_max");
			}
			printf("%5i%5i%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g\n",pl,runsCount,len,i0,s,checkSol.back(),\
				checkSolMax.back(),checkDelta.back(),checkSm.back(),checkKgDxMax.back(),checkICMax.back(),checkCCMax.back());
		}
	
	}
	
/*----------------------------------------------------------------------------------------------------------------------------
	11 - printing results
----------------------------------------------------------------------------------------------------------------------------*/	

	// eigenvalues, if required
	if (eigen) {
		mat dds_wlm = dds.block(0,0,dim*N,dim*N); // dds without Lagrange multipliers
		number eigenTol = 1.0e-16*dim*N;
		uint negEigs = 0;
		cout << "calculating eigendecomposition of dds..." << endl;
		Eigen::SelfAdjointEigenSolver<mat> eigensolver(dds_wlm);
		if (eigensolver.info()!=Eigen::Success) abort();
		//cout << "first " << 2*dim << " eigenvalues are: " << endl;
		for (uint j=0; j<8; j++) {
			if ((eigensolver.eigenvalues())[j]<-eigenTol)
				negEigs++;
			//cout << (eigensolver.eigenvalues())[j] << endl;
		}
		cout << negEigs << " negative eigenvalues found, less than " << -eigenTol << endl;
		string eigenFile = "data/nr/eigenvalues/dim_"+nts(dim)+"/K_"+nts(p.K)+"/eigenvalues_G_"+nts(p.G)+"_B_"\
							+nts(p.B)+"_M_"+nts(M)+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+".dat";
		saveVectorBinary(eigenFile,eigensolver.eigenvalues());
		printf("%12s%50s\n","eigenvalues:",((string)eigenFile).c_str());
	}
	
	// curvature, if required
	if (curvature) {
		Filename file = "data/temp/xCurvature_K_"+nts(p.K)+"_G_"+nts(p.G)+"_B_"+nts(p.B)+"_M_"+nts(M)\
						+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+".dat";
		printAsLoop(file,dim,x,N*dim);
		saveVectorAsciiAppend(file,sc_vec);
		saveVectorAsciiAppend(file,kg_vec);
		printf("%12s%50s\n","curvature:",((string)file).c_str());
	}

	//stopping clock
	time = clock() - time;
	number realtime = time/1000000.0;
	
	// results to compare to
	number s_cf = PI/p.G/p.B;
	if (!weak) s_cf -= p.G*p.G/4.0;
	else if (abs(M)>MIN_NUMBER && abs(M)<2.0) s_cf -= (2.0/p.G/p.B)*( asin(M/2.0) + (M/2.0)*sqrt(1.0-pow(M/2.0,2)) );
		
	// printing results to terminal
	printf("\n");
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%12s%14s%14s%14s%14s\n","runs","time","K","G","B","Ng","a","mu","E","len",\
		"i0","vr","s");
	printf("%8i%8.3g%8i%8.4g%8.4g%8i%8.4g%8.4g%12.5g%14.5g%14.5g%14.5g%14.5g\n",\
		runsCount,realtime,p.K,p.G,p.B,p.Ng,p.Epsi,p.Mu,M,len,i0,vr,s);
	printf("\n");
	
	
	if (checkDelta.good() && checkSol.good() && checkSolMax.good()) {
	
		// printing results to file	
		string resFile = "results/nr/nrmain_cosmos_6.dat";
		FILE* ros;
		ros = fopen(resFile.c_str(),"a");
		fprintf(ros,"%24s%24i%24i%24g%24g%24i%24g%24g%24g%24g%24g%24g%24g%24g%24g%24g%24g%24g\n",\
					timenumber.c_str(),pl,p.K,p.G,p.B,p.Ng,p.Epsi,p.Mu,M,s,(gamma0+gamma1)/2.0,\
					checkSol.back(),checkDX.back(),checkICMax.back(),checkICAvg.back(),checkKgAMax.back(),checkKgAAvg.back()\
					,checkCCMax.back());
		fclose(ros);
		printf("%12s%50s\n","results:",resFile.c_str());
		
		// printing loop to file
		Filename loopRes = filenameLoopNR<dim>(p);
		if (p.Ng>0)
			(loopRes.Extras).push_back(StringPair("Ng",nts(p.Ng)));
		saveVectorBinary(loopRes,x);
		printf("%12s%50s\n","x:",((string)loopRes).c_str());
		
	}	

	// printing extras to ascii files
	if (po!=PrintOptions::none) {
		Filename file = "data/temp/"+timenumber+"x_K_"+nts(p.K)+"_G_"+nts(p.G)+"_B_"+nts(p.B)+"_M_"+nts(M)\
							+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+"_run_"+nts(runsCount)+".dat";
		if (p.Ng>0)
				(file.Extras).push_back(StringPair("Ng",nts(p.Ng)));
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
