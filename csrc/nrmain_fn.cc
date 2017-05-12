/*
	main program to run serial n-r method to solve classical monopole worldline equations (as a function)
*/

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>
#include "analysis.h"
#include "check.h"
#include "folder.h"
#include "genloop.h"
#include "nrloop.h"
#include "print.h"
#include "simple.h"
#include "eigen_extras.h"
#include "stepper.h"

using namespace std;

//dimension
#define dim 2

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
	0 - enums
----------------------------------------------------------------------------------------------------------------------------*/

struct PrintOptions {
	enum Option { none, all, x, mds, dds, delta, curvature};
};

struct PotentialOptions {
	enum Option { original, link, exponential, dimreg, thermal, thermalDisjoint, thermalDisjointLR, external, externalDisjoint, nonrelDisjoint};
};
uint NumberPotentialOptions = 10;

struct KineticOptions {
	enum Option { saddle, s0, len};
};

struct StepperArgv {
	enum Option { none, action, entropy};
};


int nrmain_fn(int argc, vector<string> argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1 - argv, parameters etc
----------------------------------------------------------------------------------------------------------------------------*/

// argv options
bool verbose = true;
bool guess = false;
bool straight = false;
bool nonrelativistic = false;
bool step = true;
bool weak = false;
bool eigen = false;
bool includeLagrangeMultipliers = false;
bool curvature = false;
bool conservation = false;
bool old = true;
bool gaussian = false;
bool gaussianLR = false;
bool disjoint = false;
bool fixdz = false;
bool fixtlr = false;
bool fixall = false;
bool fixodt = false;
bool fixdislr = false;
bool extended = false;
bool mu_a = false;
bool auto_a = false;
bool pass = false;
bool sometests = false;
bool alltests = false; // doing alltests
bool redo = true;
bool redoErrors = true;
bool thermal2 = true;
bool onlyselfreg = false;
string baseFolder = "";
string printOpts = "";
string potOpts = "";
string kinOpts = "";
string xIn = "";
string inputsFile = "inputs4";
string stepperArgv = "";
string stepperInputsFile = "step0";
Filename stepperOutputsFile = "results/nr/stepper/step.csv";
uint steps = 1;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("guess")==0) guess = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("straight")==0) straight = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("nonrelativistic")==0 || id.compare("nonrel")==0) nonrelativistic = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("step")==0) step = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("weak")==0) weak = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("eigen")==0) eigen = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("includeLagrangeMultipliers")==0) includeLagrangeMultipliers = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("curvature")==0) curvature = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("conservation")==0) conservation = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("old")==0) old = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("gaussian")==0 || id.compare("repulsion")==0) gaussian = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("gaussianLR")==0) gaussianLR = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("mu_a")==0) mu_a = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("auto_a")==0) auto_a = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("disjoint")==0) disjoint = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("extended")==0) extended = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("pass")==0) pass = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixdz")==0) fixdz = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixtlr")==0) fixtlr = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixall")==0) fixall = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixodt")==0) fixodt = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("fixdislr")==0) fixdislr = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("base")==0) baseFolder = (string)argv[2*j+2];
		else if (id.compare("inputs")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("print")==0) printOpts = (string)argv[2*j+2];
		else if (id.compare("pot")==0 || id.compare("potential")==0) potOpts = (string)argv[2*j+2];
		else if (id.compare("kin")==0 || id.compare("kinetic")==0) kinOpts = (string)argv[2*j+2];
		else if (id.compare("stepper")==0) stepperArgv = (string)argv[2*j+2];
		else if (id.compare("steps")==0) steps = stn<uint>(argv[2*j+2]);
		else if (id.compare("stepperInputs")==0 || id.compare("stepInputs")==0) stepperInputsFile = (string)argv[2*j+2];
		else if (id.compare("stepperOutputs")==0 || id.compare("stepResults")==0) stepperOutputsFile = (string)argv[2*j+2];
		else if (id.compare("xIn")==0) xIn = (string)argv[2*j+2];
		else if (id.compare("sometests")==0 || id.compare("someTests")==0) sometests = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("alltests")==0 || id.compare("allTests")==0) alltests = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("redo")==0) redo = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("redoErrors")==0) redoErrors = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("thermal2")==0) thermal2 = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("onlyselfreg")==0) onlyselfreg = (stn<uint>(argv[2*j+2])!=0);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}
else {
	cerr << "must provide an even number of arguments after ./nrmain" << endl;
	// printing argv
	for (int j=0; j<argc; j++) {
		cerr << argv[j];
		if (j<(argc-1))
			cerr << " ";
		else
			cerr << endl;
	}
	return 1;
}

if (verbose) {
	// printing argv
	for (int j=0; j<argc; j++) {
		cout << argv[j];
		if (j<(argc-1))
			cout << " ";
		else
			cout << endl;
	}
}

PrintOptions::Option po = PrintOptions::none;
if (!printOpts.empty()) {
	if (printOpts.compare("all")==0) {
		po = PrintOptions::all;
		curvature = true;
	}
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
	else if (printOpts.compare("none")==0) {
		po = PrintOptions::none;
	}
	else
		cerr << "print options not understood: " << printOpts << endl;
}

PotentialOptions::Option poto = PotentialOptions::original;
if (!potOpts.empty()) {
	if (potOpts.compare("original")==0)
		poto = PotentialOptions::original;
	else if (potOpts.compare("link")==0)
		poto = PotentialOptions::link;
	else if (potOpts.compare("exponential")==0 || potOpts.compare("exp")==0)
		poto = PotentialOptions::exponential;
	else if (potOpts.compare("dimreg")==0)
		poto = PotentialOptions::dimreg;
	else if (potOpts.compare("thermal")==0)
		poto = PotentialOptions::thermal;
	else if (potOpts.compare("thermalDisjoint")==0) {
		poto = PotentialOptions::thermalDisjoint;
		disjoint = true;
	}
	else if (potOpts.compare("thermalDisjointLR")==0) {
		poto = PotentialOptions::thermalDisjointLR;
		disjoint = true;
	}
	else if (potOpts.compare("external")==0)
		poto = PotentialOptions::external;
	else if (potOpts.compare("externalDisjoint")==0) {
		poto = PotentialOptions::externalDisjoint;
		disjoint = true;	
	}
	else if (potOpts.compare("nonrelDisjoint")==0) {
		poto = PotentialOptions::nonrelDisjoint;
		disjoint = true;	
	}
	else {
		cerr << "potential options not understood: " << potOpts << endl;
		return 1;
	}
}
if (!disjoint)
	gaussianLR = false;
else
	onlyselfreg = false;
if (gaussianLR)
	gaussian = true;
if (onlyselfreg)
	thermal2 = false;
if (thermal2)
	gaussian = true;

if (disjoint) { // some defaults to save time typing, not necessary, feel free to remove this
	fixdz = true;
	fixtlr = true;
	fixdislr = true;
}

StringPair potExtras("pot",nts((int)poto));
if ((int)poto<5 && !thermal2 && !onlyselfreg)
	potExtras.second = nts(2*(int)poto + (int)gaussian);
else if (!thermal2 && !onlyselfreg)
	potExtras.second = nts(10 + 3*((int)poto-5) + (int)gaussian + (int)gaussianLR);
else if (!disjoint && !onlyselfreg)
	potExtras.second = nts(13);
else if (!disjoint && onlyselfreg)
	potExtras.second = nts(16);
else if (poto == PotentialOptions::thermalDisjoint) {
	if (gaussian)
		potExtras.second = nts(14); // this is all a bit messy now but whatever
	else
		potExtras.second = nts(17);
}
else if (poto == PotentialOptions::thermalDisjointLR) {
	if (gaussianLR)
		potExtras.second = nts(15);
	else 
		potExtras.second = nts(18);
}
		
	
KineticOptions::Option kino = KineticOptions::saddle;
StringPair kinExtras("kin",nts((int)kino));
if (!kinOpts.empty()) {
	if (kinOpts.compare("saddle")==0)
		kino = KineticOptions::saddle;
	else if (kinOpts.compare("s0")==0)
		kino = KineticOptions::s0;
	else if (kinOpts.compare("len")==0 || potOpts.compare("length")==0)
		kino = KineticOptions::len;
	else {
		cerr << "potential options not understood/available: " << kinOpts << endl;
		return 1;
	}
	kinExtras.second = nts((int)kino);
}

// stepper
StepperArgv::Option stepargv = StepperArgv::none;
if (!stepperArgv.empty()) {
	if (stepperArgv.compare("action")==0 || stepperArgv.compare("s")==0 || stepperArgv.compare("S")==0) {
		stepargv = StepperArgv::action;
		(stepperOutputsFile.Extras).push_back(StringPair("const","action"));
		if (verbose)
			cout << "stepping with constant action" << endl;
	}
	else if (stepperArgv.compare("entropy")==0 || stepperArgv.compare("sigma")==0) {
		stepargv = StepperArgv::entropy;
		(stepperOutputsFile.Extras).push_back(StringPair("const","entropy"));
		if (verbose)
			cout << "stepping with constant entropy" << endl;
	}
	else if (stepperArgv.compare("none")==0)
		stepargv = StepperArgv::none;
	else {
		cerr << "stepper options not understood/available: " << stepperArgv << endl;
		return 1;
	}
}

if (fixall)
	fixtlr = false;

// parameters
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min, pold = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file: " << inputsFile << endl;
	return 1;
}

// initializing stepper
StepperOptions stepOpts;
stepOpts.tol = 1.0;
Point2d point;
if (stepargv!=StepperArgv::none) {
	{
		ifstream is;
		is.open(stepperInputsFile.c_str());
		if (is.good()) {
			is >> stepOpts;
			is.close();
			(stepperOutputsFile.Extras).push_back(StringPair("tol",nts(stepOpts.tol)));
			(stepperOutputsFile.Extras).push_back(StringPair("aim",nts(stepOpts.aim)));
		}
		else {
			cerr << "Error: cannot open stepper inputs file, " << stepperInputsFile << endl;
			return 1;
		}
	}
	if (poto==PotentialOptions::thermal || disjoint) {
		point(p.B,p.T);
		//stepOpts.epsi_x *= (abs(p.B)>MIN_NUMBER? p.B: 1.0);
		//stepOpts.epsi_y *= (abs(p.T)>MIN_NUMBER? p.T: 1.0);
	}
	else {
		point(p.B,p.P4);
		//stepOpts.epsi_x *= (abs(p.B)>MIN_NUMBER? p.B: 1.0);
		//stepOpts.epsi_y *= (abs(p.P4)>MIN_NUMBER? p.P4: 1.0);
	}
}
Stepper stepper(stepOpts,point);
if (stepperOutputsFile.exists() && stepargv!=StepperArgv::none) {
	stepper.load(stepperOutputsFile);
	if (poto==PotentialOptions::thermal || disjoint) {
		p.B = stepper.x();
		p.T = stepper.y();
		pold.B = (stepper.lastStep()).X;
		pold.T = (stepper.lastStep()).Y;
	}
	else {
		p.B = stepper.x();
		p.P4 = stepper.y();
		pold.B = (stepper.lastStep()).X;
		pold.P4 = (stepper.lastStep()).Y;
	}
}

// results
string resultsFile = (pass? "results/nr/nr7pass.csv":"results/nr/nr7.csv");
uint idSizeResults = 3, datumSizeResults = 17;
vector<string> idCheck(idSizeResults);
idCheck[idSizeResults-1] = potExtras.second;
NewtonRaphsonData results(resultsFile,idSizeResults,datumSizeResults);

// errors
string errorsFile = "results/nr/nr7error.csv";
uint idSizeErrors = 4, datumSizeErrors = 13;
vector<string> idCheckErrors(idSizeErrors);
idCheckErrors[idSizeErrors-1] = potExtras.second;
NewtonRaphsonData errors(errorsFile,idSizeErrors,datumSizeErrors);

// timenumber
string timenumber = currentDateTime();
cout << "timenumber: " << timenumber << ", pot: " << potExtras.second << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	2 - beginning parameter loop
----------------------------------------------------------------------------------------------------------------------------*/

// parameter loops
uint Npl = (stepargv==StepperArgv::none? pr.totalSteps(): steps);
if (verbose)
	cout << "looping over " << Npl << " steps" << endl;

// starting loop
for (uint pl=0; pl<Npl; pl++) {

	// doing things before parameter step
	Filename loadFile, stepFile;
	
	// stepping parameters
	if (pl>0) {
		if (stepargv==StepperArgv::none) {
			p = pr.position(pl);
			pold = pr.neigh(pl);
		}
		else {
			// getting step base
			if (poto==PotentialOptions::thermal || disjoint) {
				p.B = stepper.x();
				p.T = stepper.y();
				pold.B = (stepper.lastStep()).X;
				pold.T = (stepper.lastStep()).Y;
			}
			else {
				p.B = stepper.x();
				p.P4 = stepper.y();
				pold.B = (stepper.lastStep()).X;
				pold.P4 = (stepper.lastStep()).Y;
			}
		}
		if (mu_a) {
			p.Mu = p.Epsi;
			pold.Mu = pold.Epsi;
		}
	}
	
	// checking if have results already
	if (!redo && results.find(idCheck,p)) {
		if (verbose) {
			cout << "result found in " << resultsFile << " for pl = " << pl << ", ";
			cout << "continuing to next step" << endl;
		}
		continue;
	}
	if (!redoErrors && errors.find(idCheckErrors,p)) {
		if (verbose) {
			cout << "result found in " << errorsFile << " for pl = " << pl << ", ";
			cout << "continuing to next step" << endl;
		}
		continue;
	}
	
	// getting step file
	if (pl>0) {	
		if (poto==PotentialOptions::thermal || disjoint) {
			if (p.T==0) {
				cerr << "must have T!=0 for thermal runs" << endl;
				return 1;
			}
			stepFile = filenameThermalNR<dim>(pold,baseFolder);
		}
		else
			stepFile = filenameLoopNR<dim>(pold,baseFolder);
		if (weak)
			(stepFile.Extras).push_back(StringPair("weak","1"));
		if (poto!=PotentialOptions::original || gaussian)
			(stepFile.Extras).push_back(potExtras);
		if (kino!=KineticOptions::saddle)
			(stepFile.Extras).push_back(kinExtras);
	}
		
	
	// defining some derived parameters	
	uint N = pow(2,p.K);
	uint zm = dim;
	if (fixdz || fixodt)
		zm += 1+(uint)fixdislr;
	if (fixall)
		zm += dim;
	else if (fixtlr)
		zm += 1;
	uint NT = N*dim+zm;
	number R = 1.0;
	 ////////////////////////////////////// dimension dependent bit ///////////////////////////////////////////////////
	Point<dim> P;
	P[0] = p.P3;
	P[1] = p.P4;
	//P[2] = p.P3;
	//P[3] = p.P4;
	 ////////////////////////////////////// dimension dependent bit ///////////////////////////////////////////////////
	number E = p.P4;

/*----------------------------------------------------------------------------------------------------------------------------
	3 - defining quantitites
----------------------------------------------------------------------------------------------------------------------------*/

	//defining a time and starting the clock
	clock_t time;
	time = clock();
	
	// defining checks
	Check checkSol("solution",1.0e-7);
	Check checkSolMax("solution max",1.0e-6);
	Check checkSolZM("solution zero modes",1.0e-3);
	Check checkDelta("delta",1.0);
	Check checkDeltaMax("delta max",1.0);
	Check checkSm("smoothness",1.0);
	Check checkDX("dx<<a",2.0e-1);
	Check checkStraight("angle_neigh<<1",2.0e-1);
	Check checkICMax("ic_max<<1",5.0e-1);
	Check checkCCMax("cc_max<<1",5.0e-1);
	Check checkKgAMax("a*kg_max<<1",5.0e-1);
	Check checkKgDxMax("dx*kg_max<<1",5.0e-1);
	Check checkAccAMax("a*acc<<1",5.0e-1);
	Check checkSym("symmetric",1.0e-16*NT*NT);
	Check checkSymMax("symmetric max",1.0e-14);
	Check checkInv("inverse",1.0e-16*NT*NT*NT);
	Check checkNegEigenvalue("negative eigenvalue",0.1);
	Check checkNegEigenvector("negative eigenvector",0.1);
	Check checkGamma("gamma",0.1);
	Check checkJs("Js conservation",1.0e-3);
	Check checkXMirror("x mirror symmetry",1.0e-2);
	Check checkXRotation("x rotation symmetry",1.0e-16*NT*NT);
	Check checkMDSMirror("mds mirror symmetry",1.0e-2);
	Check checkMDSRotation("mds rotation symmetry",1.0e-16*NT*NT);
	Check checkDeltaMirror("delta mirror symmetry",1.0e-2);
	Check checkDeltaRotation("delta rotation symmetry",1.0e-16*NT*NT);
	Check checkEAgree("energies agree",1.0e-3);
	
	// defining scalar quantities
	number len, i0, kinetic = 0.0, s, sm, v = 0.0, vr = 0.0, erg, ergThermal = 0.0, fgamma, gamma, angle_neigh, zmax, zmin, tmax, ic_max, cc_max, kg_max, acc_max;
	
	// defining vector and matrix quantities
	vec x(N*dim);
	vec mds(NT);
	vec delta(NT);
	mat dds(NT,NT);
	
	// curvature
	vec sc_vec;
	vec kg_vec;
	vec acc_vec;
	
	// conserved quantities
	vec Js(N);
	
	// defining xLoop
	Loop<dim> xLoop(p.K,0);
	
	// x file
	if (pl==0 || !step || pass) {
		if (!xIn.empty()) {
			loadFile = xIn;
		}
		else if (guess) {
			if (!disjoint) {
				loadFile = baseFolder+"data/lemon/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R)+"_E_"+nts(E)+"_rank_0.dat";
				if (!loadFile.exists())
					loadFile = baseFolder+"data/circle/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R)+"_rank_0.dat";
			}
			else {
				if (nonrelativistic)
					loadFile = baseFolder+"data/highTemp/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/highTemp_kappa_"+nts(pow(p.G,3)*p.B)\
					+"_T_"+nts(p.T)+"_rank_0.dat";
				else 
					loadFile = baseFolder+"data/cosDisjoint/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)\
						+"_T_"+nts(p.T)+"_mu_"+nts(p.Mu)+"_lambda_"+nts(p.Lambda)+"_rank_0.dat";
				if (!loadFile.exists() || straight) {
					loadFile = baseFolder+"data/straightDisjoint/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)\
					+"_T_"+nts(p.T)+"_rank_0.dat";
					if (extended)
						(loadFile.Extras).push_back(StringPair("Lambda",nts(p.Lambda)));
				}
			}
		}
		else {
			if (poto==PotentialOptions::thermal || disjoint)
				loadFile = filenameThermalNR<dim>(p,baseFolder);
			else
				loadFile = filenameLoopNR<dim>(p,baseFolder);
			if (weak)
				(loadFile.Extras).push_back(StringPair("weak","1"));
			if (poto!=PotentialOptions::original || gaussian)
				(loadFile.Extras).push_back(potExtras);
			if (kino!=KineticOptions::saddle)
				(loadFile.Extras).push_back(kinExtras);
			if (!loadFile.exists() && (poto!=PotentialOptions::original || gaussian)) {
				loadFile = filenameLoopNR<dim>(p,baseFolder);
				(loadFile.Extras).push_back(potExtras);
			}
			if (!loadFile.exists() && poto==PotentialOptions::thermal && !thermal2) {
				loadFile = filenameThermalNR<dim>(p,baseFolder);
				StringPair potExtrasAlt("pot","13");
				(loadFile.Extras).push_back(potExtrasAlt);
				if (!loadFile.exists()) {
					loadFile = filenameThermalNR<dim>(p,baseFolder);
					StringPair potExtrasAlt("pot","8");
					(loadFile.Extras).push_back(potExtrasAlt);
				}
			}
			else if (!loadFile.exists() && poto==PotentialOptions::thermalDisjoint && !thermal2) {
				loadFile = filenameThermalNR<dim>(p,baseFolder);
				StringPair potExtrasAlt("pot","");
				int offset = (gaussianLR? -2: 2);
				potExtrasAlt.second = nts((int)(stn<int>(potExtras.second)+offset));
				(loadFile.Extras).push_back(potExtrasAlt);
			}
			else if (!loadFile.exists() && poto==PotentialOptions::thermal && thermal2) {
				loadFile = filenameThermalNR<dim>(p,baseFolder);
				StringPair potExtrasAlt("pot","");
				potExtrasAlt.second = nts(9);
				(loadFile.Extras).push_back(potExtrasAlt);
				if (!loadFile.exists()) {
					loadFile = filenameThermalNR<dim>(p,baseFolder);
					StringPair potExtrasAlt("pot","");
					potExtrasAlt.second = nts(8);
					(loadFile.Extras).push_back(potExtrasAlt);
				}
			}
			else if (!loadFile.exists() && poto==PotentialOptions::thermalDisjoint && thermal2) {
				loadFile = filenameThermalNR<dim>(p,baseFolder);
				StringPair potExtrasAlt("pot","");
				potExtrasAlt.second = nts(12);
				(loadFile.Extras).push_back(potExtrasAlt);
				if (!loadFile.exists()) {
					loadFile = filenameThermalNR<dim>(p,baseFolder);
					StringPair potExtrasAlt("pot","");
					potExtrasAlt.second = nts(10);
					(loadFile.Extras).push_back(potExtrasAlt);
				}
			}
			else if (!loadFile.exists() && poto==PotentialOptions::thermalDisjointLR  && thermal2) { // doesn't work in many circumstances, just quick and dirty
				loadFile = filenameThermalNR<dim>(p,baseFolder);
				StringPair potExtrasAlt("pot","");
				potExtrasAlt.second = nts(14);
				(loadFile.Extras).push_back(potExtrasAlt);
				if (!loadFile.exists()) {
					loadFile = filenameThermalNR<dim>(p,baseFolder);
					StringPair potExtrasAlt("pot","");
					potExtrasAlt.second = nts(12);
					(loadFile.Extras).push_back(potExtrasAlt);
					if (!loadFile.exists()) {
						loadFile = filenameThermalNR<dim>(p,baseFolder);
						StringPair potExtrasAlt("pot","");
						potExtrasAlt.second = nts(10);
						(loadFile.Extras).push_back(potExtrasAlt);
					}
				}
			}
			else if (!loadFile.exists() && !thermal2) {
				loadFile = filenameLoopNR<dim>(p,baseFolder);
				StringPair potExtrasAlt("pot","");
				int offset = (gaussian? -1:+1);
				potExtrasAlt.second = nts((int)(stn<int>(potExtras.second)+offset));
				(loadFile.Extras).push_back(potExtrasAlt);
			}
			if (!loadFile.exists() && (poto==PotentialOptions::thermal || disjoint))
				loadFile = filenameThermalNR<dim>(p,baseFolder);
			else if (!loadFile.exists())
				loadFile = filenameLoopNR<dim>(p,baseFolder);
			if (!loadFile.exists()) {
				if (pl>0)
					loadFile = stepFile;
					if (!loadFile.exists() && (poto==PotentialOptions::thermal || disjoint))
						loadFile = filenameLoopNR<dim>(pold,baseFolder);
					else if (!loadFile.exists())
						loadFile = filenameThermalNR<dim>(pold,baseFolder);
				if (!loadFile.exists() && old)
					loadFile = baseFolder+"data/nr/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(p.P4)\
		+"_a_"+nts(p.Epsi)+".dat";
		
				if (!loadFile.exists()) {
					if (!disjoint) {
					loadFile = baseFolder+"data/lemon/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R)+"_E_"+nts(E)+"_rank_0.dat";
						if (!loadFile.exists())
							loadFile = baseFolder+"data/circle/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_R_"+nts(R)+"_rank_0.dat";
					}
					else {
						if (nonrelativistic)
							loadFile = baseFolder+"data/highTemp/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/highTemp_kappa_"+nts(pow(p.G,3)*p.B)\
							+"_T_"+nts(p.T)+"_rank_0.dat";
						else 
							loadFile = baseFolder+"data/cosDisjoint/loops/dim_"+nts(dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)\
								+"_T_"+nts(p.T)+"_mu_"+nts(p.Mu)+"_lambda_"+nts(p.Lambda)+"_rank_0.dat";
						if (!loadFile.exists() || straight) {
							loadFile = baseFolder+"data/straightDisjoint/loops/dim_"+nts(dim)+"/K_"+nts(p.K)\
							+"/loop_kappa_"+nts(pow(p.G,3)*p.B)+"_T_"+nts(p.T)+"_rank_0.dat";
							if (extended)
								(loadFile.Extras).push_back(StringPair("Lambda",nts(p.Lambda)));
						}
					}
				}
			}
		}
			
	}
	else {
		loadFile = stepFile;
		guess = false;
	}
	// check if file exists
	if (!loadFile.exists()) {
		cerr << "nrmain error: " << loadFile << " doesn't exist on pl = " << pl << ", moving to next parameter loop" << endl;
		continue; ///////// CONTINUE STATEMENT IF FILE DOESN'T EXIST
	}
	if (verbose) {
		cout << "loading loops from:" << endl;
		cout << loadFile << endl;
	}
	
	// loading x
	loadVectorBinary(loadFile,x);
	if (pold.K!=p.K && loadFile==stepFile) {
		Loop<dim> interpOld(pold.K,0), interpNew(p.K,0);
		vectorToLoop(x,interpOld);
		interpolate(interpOld,interpNew);
		loopToVector(interpNew,x);
	}
	if (x.size()<NT || pl==0) {
		x.conservativeResize(NT);
		for (uint mu=0; mu<zm; mu++)
			x[N*dim+mu] = 1.0e-4;
	}
	else if (x.size()>NT)
		x.conservativeResize(NT);
		
	if (disjoint) {
		// inverting lhs for disjoint topology, if necessary
		if (x[(N-1)*dim+dim-1]>x[(N-2)*dim+dim-1]) { // switched again
			number tempx;
			for (uint j=0; j<N/4; j++) {
				for (uint mu=0; mu<dim; mu++) {
					tempx = x[(N/2+j)*dim+mu];
					x[(N/2+j)*dim+mu] = x[(N-1-j)*dim+mu];
					x[(N-1-j)*dim+mu] = tempx;
				}
			}
		}
		// stretching if changing T and disjoint
		if (abs(p.T-pold.T)>MIN_NUMBER && pl>0 && step) {
			for (uint j=0; j<N; j++) {
				x[dim*j+dim-1] *= pold.T/p.T;
			}
		}
	}
	
	// if auto_a, checking "a" is small enough
	if (auto_a) {
		vectorToLoop(x,xLoop);
		len = L(xLoop);
		number dx = len/(number)N;
		if ((dx/p.Epsi)>p.Lambda && p.Lambda>MIN_NUMBER && dx>MIN_NUMBER) {
			number a_new = sigFig(dx/p.Lambda,2.0);
			if (abs(a_new-p.Epsi)>MIN_NUMBER) {
				if (verbose)
					cout << "changed a: a_old = " << p.Epsi << ", a_new = " << a_new << ", using Lambda = " << p.Lambda << endl;
				(pr.Min).Epsi 	= a_new;
				(pr.Max).Epsi 	= a_new;
				p.Epsi 		= a_new;
				pold.Epsi 	= a_new;
			}
		}
	}
	
	//defining some quantities used to stop n-r loop
	uint runsCount = 0;
	uint minRuns = 1;
	uint maxRuns = 100;
	
/*----------------------------------------------------------------------------------------------------------------------------
	4 - beginning newton-raphson loop
----------------------------------------------------------------------------------------------------------------------------*/

	bool passThrough = false;
	while ((!checkSol.good() || !checkSolMax.good() || runsCount<minRuns) && runsCount<maxRuns && !passThrough)  {
		runsCount++;
		passThrough = pass;
		
		// initializing to zero
		mds = Eigen::VectorXd::Zero(NT);
		dds = Eigen::MatrixXd::Zero(NT,NT);
		Js = Eigen::VectorXd::Zero(N);
		len = 0.0, kinetic = 0.0, i0 = 0.0, v = 0.0, fgamma = 0.0, gamma = 0.0, angle_neigh = 0.0, zmax = 0.0, tmax = 0.0, erg = 0.0;
		zmin = 1.0e16;
		ic_max = 0.0, cc_max = 0.0, kg_max = 0.0, acc_max = 0.0;
		if (curvature || alltests) {
			sc_vec = Eigen::VectorXd::Zero(N);
			kg_vec = Eigen::VectorXd::Zero(N);
			acc_vec = Eigen::VectorXd::Zero(N);
		}
				
		// loading x to xLoop - messier than it should be (should work with either a vec or a Loop really)
		vectorToLoop(x,xLoop);

/*----------------------------------------------------------------------------------------------------------------------------
	5 - assigning mdx and ddx
----------------------------------------------------------------------------------------------------------------------------*/
		
		// scalar coefficients
		uint j, k, mu, nu;
		number mgb = -1.0; // not -p.G*p.B as scaled loops
		number g, dm, cusp_scale;
		number dim_reg_scale = 0.0, d_dim_reg = 0.0;
		number n = 0.0;
		number repulsion_scale, repulsion = 0.0;
		number ic_scale = 1.0;
		number cc_scale = 1.0;
		number kg_scale = 1.0;
		number Js_scale = 4.0*N;
		number s0, sqrt4s0; 
		number s0_scale = (abs(p.T)>MIN_NUMBER? 1.0/p.T: 1.0);
		number beta = ((p.T)>sqrt(MIN_NUMBER)? 1.0/(p.T): 1.0/sqrt(MIN_NUMBER)); // this is 1/eta
		len = (disjoint? LDisjoint(xLoop,beta): L(xLoop));
		number acc_scale = 1.0/pow(len,2);
		number sigma = 1.0;
		if (poto==PotentialOptions::original) {
			s0 = S0(xLoop);
			sqrt4s0 = 2.0*sqrt(s0);
			g = pow(p.G,3)*p.B/8.0/PI/PI;
			dm = -g*PI/p.Epsi;
			cusp_scale = -g*2.0*log(p.Mu/p.Epsi);
			repulsion_scale = -g*sqrt(PI)/p.Epsi/p.Epsi;
		}
		else if (poto==PotentialOptions::link) {
			s0 = S0(xLoop);
			sqrt4s0 = 2.0*sqrt(s0);
			g = pow(p.G,3)*p.B/8.0/PI/PI;
			dm = -g*PI/p.Epsi;
			cusp_scale = -g*2.0*log(p.Mu/p.Epsi);
			repulsion_scale = -g*sqrt(PI)/p.Epsi/p.Epsi;
		}
		else if (poto==PotentialOptions::exponential) {
			s0 = S0(xLoop);
			sqrt4s0 = 2.0*sqrt(s0);
			g = pow(p.G,3)*p.B/8.0/PI/PI;
			dm = -g*sqrt(PI)/p.Epsi;
			cusp_scale = -g*2.0*log(p.Mu/p.Epsi);
			repulsion_scale = -g*sqrt(PI)/p.Epsi/p.Epsi; // just added *sqrt(PI) to this line
		}
		else if (poto==PotentialOptions::dimreg) {
			s0 = S0(xLoop);
			sqrt4s0 = 2.0*sqrt(s0);
			g = pow(p.G,3)*p.B*(pow(p.Lambda,-p.Epsi)*pow(PI,-2.0+p.Epsi/2.0)*gsl_sf_gamma(2.0-p.Epsi/2.0)/4.0/(2.0-p.Epsi));
			dm = 0.0;
			cusp_scale = -g*2.0/p.Epsi; //not sure about the factor of 2.0 here
			repulsion_scale = 0.0;
			dim_reg_scale = -g*2.0*gsl_sf_zeta(2.0-p.Epsi);
		}
		else if (poto==PotentialOptions::thermal || poto==PotentialOptions::thermalDisjoint || poto==PotentialOptions::thermalDisjointLR) {
			g = pow(p.G,3)*p.B/8.0/PI/PI;
			dm = -g*PI/p.Epsi;
			cusp_scale = -g*2.0*log(p.Mu/p.Epsi);
			repulsion_scale = -g*sqrt(PI)/p.Epsi/p.Epsi;
			if (disjoint)
				s0 = S0Disjoint(xLoop,beta);
			else
				s0 = S0(xLoop);
			sqrt4s0 = 2.0*sqrt(s0);
		}
		else if (poto==PotentialOptions::external || poto==PotentialOptions::externalDisjoint) {
			n = -1.0;
			g = -pow(p.G,3)*p.B/4.0/PI/4.0;
			dm = 0.0;
			cusp_scale = 0.0;
			repulsion_scale = 0.0;
			if (disjoint)
				s0 = S0Disjoint(xLoop,beta);
			else
				s0 = S0(xLoop);
			sqrt4s0 = 2.0*sqrt(s0);
		}
		else if (poto==PotentialOptions::nonrelDisjoint) {
			n = -1.0;
			g = -pow(p.G,3)*p.B/4.0/PI;
			dm = 0.0;
			cusp_scale = 0.0;
			repulsion_scale = 0.0;
			if (disjoint)
				s0 = S0Disjoint(xLoop,beta);
			else
				s0 = S0(xLoop);
			sqrt4s0 = 2.0*sqrt(s0);
		}
		
		Point<dim> P0;
		
		// vector carrying nonlocal remainder in noether charge
		vec PRmu = Eigen::VectorXd::Zero(N*dim);
		
		// bulk
		for (j=0; j<N; j++) {
			
			//S0(j, xLoop, s0norm, s0norm);
			if (!disjoint) {
				//L		(j, xLoop, 1.0, len);
				I0		(j, xLoop, mgb, i0);
				S0		(j, xLoop, Js_scale, Js(j));
				if (poto==PotentialOptions::external)
					In(j, xLoop, n, g, v);
			}
			else {
				//LDisjoint (j, xLoop, beta, 1.0, len);
				I0Disjoint(j, xLoop, beta, mgb, i0);
				S0Disjoint(j, xLoop, beta, Js_scale, Js(j));
				if (poto==PotentialOptions::externalDisjoint)
					InDisjoint(j, xLoop, n, beta, g, v);
			}
			
			if (poto==PotentialOptions::dimreg)
				DistPow	(j, xLoop, p.Epsi, dim_reg_scale, d_dim_reg);
			
			//curvatures
			if (curvature || alltests) {
				if (!disjoint) {
					InlineCurvatureMax(j, xLoop, ic_scale, sc_vec[j]);
					KGMaxPlane(j, xLoop, kg_scale, kg_vec[j]);
					AccMax(j, xLoop, acc_scale, acc_vec[j]);
				}
				else {
					InlineCurvatureMaxDisjoint(j, xLoop, beta, ic_scale, sc_vec[j]);
					KGMaxPlaneDisjoint(j, xLoop, beta, kg_scale, kg_vec[j]);
					AccMaxDisjoint(j, xLoop, beta, acc_scale, acc_vec[j]);
				}
			}
			if (!(P^=P0)) {
				if (!disjoint) {
					InlineCurvatureMax(j, xLoop, 0, N/2-1, ic_scale, ic_max);
					CuspCurvatureMax(j, xLoop, 0, N/2-1, cc_scale, ic_max);
					KGMaxPlane(j, xLoop, 0, N/2-1, kg_scale, kg_max);
					AccMax(j, xLoop, 0, N/2-1, acc_scale, acc_max);
				}
				else {
					InlineCurvatureMaxDisjoint(j, xLoop, 0, N/2-1, beta, ic_scale, ic_max);
					CuspCurvatureMaxDisjoint(j, xLoop, 0, N/2-1, beta, cc_scale, ic_max);
					KGMaxPlaneDisjoint(j, xLoop, 0, N/2-1, beta, kg_scale, kg_max);
					AccMaxDisjoint(j, xLoop, 0, N/2-1, beta, acc_scale, acc_max);
				}
			} 
			else {
				if (!disjoint) {
					InlineCurvatureMax(j, xLoop, ic_scale, ic_max);
					CuspCurvatureMax(j,xLoop,cc_scale,ic_max);
					KGMaxPlane(j, xLoop, kg_scale, kg_max);
					AccMax(j, xLoop, acc_scale, acc_max);
				}
				else {
					InlineCurvatureMaxDisjoint(j, xLoop, beta, ic_scale, ic_max);
					CuspCurvatureMaxDisjoint(j,xLoop, beta,cc_scale,ic_max);
					KGMaxPlaneDisjoint(j, xLoop, beta, kg_scale, kg_max);
					AccMaxDisjoint(j, xLoop, beta, acc_scale, acc_max);
				}
			}
		
			for (mu=0; mu<dim; mu++) {
			
				// free particle
				if (!disjoint) {
					if (kino==KineticOptions::saddle) {
						mdsqrtS0_nr(j,mu,xLoop,sqrt4s0,1.0,mds);
						if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgsqrtS0_nr(xLoop, j, mu, sqrt4s0, sigma*1.0, erg);
						}
					}
					else if (kino==KineticOptions::s0) {
						mdS0_nr(j,mu,xLoop,s0_scale,mds);
						if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgS0_nr(xLoop, j, mu, sigma*s0_scale, erg);
						}
					}
					else if (kino==KineticOptions::len) {
						mdL_nr(j,mu,xLoop,1.0,mds);
						if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgL_nr(xLoop, j, mu, sigma*1.0, erg);
						}
					}
				}
				else {
					if (kino==KineticOptions::saddle) {
						mdsqrtS0Disjoint_nr(j,mu,xLoop,sqrt4s0,beta,1.0,mds);
						if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgsqrtS0Disjoint_nr(xLoop, j, mu, sqrt4s0, beta, sigma*1.0, erg);
						}
					}
					else if (kino==KineticOptions::s0) {
						mdS0Disjoint_nr(j,mu,xLoop,beta,s0_scale,mds);
						if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgS0Disjoint_nr(xLoop, j, mu, beta, sigma*s0_scale, erg);
						}
					}
					else if (kino==KineticOptions::len) {
						mdLDisjoint_nr(j,mu,xLoop,beta,1.0,mds);
						if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgLDisjoint_nr(xLoop, j, mu, beta, sigma*1.0, erg);
						}
					}
				}
				
				// external field
				if (!disjoint) {
					mdI0_nr(j, mu, xLoop, mgb, mds);
					if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
						sigma = sign((xLoop[j])[dim-2]);
						ErgI0_nr(xLoop, j, mu, sigma*mgb, erg);
					}
					if (poto==PotentialOptions::external) {
						mdIn_nr(j, mu, xLoop, n, g, mds);
						if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgIn_nr(xLoop, j, mu, n, sigma*g, erg);
						}
					}
				}
				else {
					mdI0Disjoint_nr(j,mu,xLoop,beta,mgb,mds);
					if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
						sigma = sign((xLoop[j])[dim-2]);
						ErgI0Disjoint_nr(xLoop, j, mu, beta, sigma*mgb, erg);
					}
					if (poto==PotentialOptions::externalDisjoint) {
						mdInDisjoint_nr(j, mu, xLoop, n, beta, g, mds);
						if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgInDisjoint_nr(xLoop, j, mu, n, beta, sigma*g, erg);
						}
					}
					else if (poto==PotentialOptions::nonrelDisjoint) {
						VnonrelrDisjoint(j, xLoop, beta, p.Epsi, g, v); //
						mdVnonrelrDisjoint_nr(j, mu, xLoop, beta, p.Epsi, g, mds);
					}
				}
				
				if (!weak && !gaussian) {
					// self-energy regularisation
					if (poto!=PotentialOptions::dimreg && !disjoint) {
						mdL_nr(j,mu,xLoop,dm,mds);
						if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgL_nr(xLoop, j, mu, sigma*dm, erg);
						}
					}
					else if (disjoint) {
						mdLDisjoint_nr(j,mu,xLoop,beta,dm,mds);
						if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
							sigma = sign((xLoop[j])[dim-2]);
							ErgLDisjoint_nr(xLoop, j, mu, beta, sigma*dm, erg);
						}
					}
					else {
						mdDistPow_nr(j, mu, xLoop, p.Epsi, dim_reg_scale, mds);
					}
				}
				
				for (k=0; k<N; k++) {
				
					if (mu==0) {
						if (poto==PotentialOptions::original)
							Vor(j, k, xLoop, p.Epsi, g, v);
						else if (poto==PotentialOptions::link)
							Vlr(j, k, xLoop, p.Epsi, g, v);
						else if (poto==PotentialOptions::exponential)
							Ver(j, k, xLoop, p.Epsi, g, v);
						else if (poto==PotentialOptions::dimreg)
							Vdr(j, k, xLoop, p.Epsi, g, v);
						else if (poto==PotentialOptions::thermal) {
							if (onlyselfreg)
								Vth2r(j, k, xLoop, beta, p.Epsi, g, v);
							else
								Vthr(j, k, xLoop, beta, p.Epsi, g, v);
						}
						else if (poto==PotentialOptions::thermalDisjoint)
							VthrDisjoint(j, k, xLoop, beta, p.Epsi, g, v);
						else if (poto==PotentialOptions::thermalDisjointLR)
							VthrDisjointLR(j, k, xLoop, beta, p.Epsi, g, v);
							
						if (gaussian) {
							if (!disjoint && poto!=PotentialOptions::thermal)
								Gaussian(j, k, xLoop, p.Epsi, repulsion_scale, repulsion);
							else if (disjoint && poto!=PotentialOptions::thermalDisjoint && poto!=PotentialOptions::thermalDisjointLR) {
								if (gaussianLR)
									GaussianLRDisjoint(j, k, xLoop, beta, p.Epsi, repulsion_scale, repulsion);
								else
									GaussianDisjoint(j, k, xLoop, beta, p.Epsi, repulsion_scale, repulsion);
							}
							else if (poto==PotentialOptions::thermal) {
								if (thermal2)
									GaussianThermal2(j, k, xLoop, beta, p.Epsi, repulsion_scale, repulsion);
								else
									Gaussian(j, k, xLoop, p.Epsi, repulsion_scale, repulsion);
							}
							else if (poto==PotentialOptions::thermalDisjoint || poto==PotentialOptions::thermalDisjointLR) {
								if (gaussianLR) {
									if (thermal2)
										GaussianThermalLRDisjoint(j, k, xLoop, beta, p.Epsi, repulsion_scale, repulsion);
									else
										GaussianThermalLRDisjoint(j, k, xLoop, beta, p.Epsi, repulsion_scale, repulsion);
								}
								else {
									if (thermal2)
										GaussianThermal2Disjoint(j, k, xLoop, beta, p.Epsi, repulsion_scale, repulsion);
									else
										GaussianThermalDisjoint(j, k, xLoop, beta, p.Epsi, repulsion_scale, repulsion);
								}
							}
						}
							
						MaxXn(j, k, xLoop, dim-2, 1.0, zmax);
						MinXnDisjoint(j, k, xLoop, dim-2, 1.0, zmin);
						MaxXn(j, k, xLoop, dim-1, 1.0, tmax);
					}
					
					// dynamical field
					if (!weak) {
						if (poto==PotentialOptions::original) {
							mdVor_nr(j, mu, k, xLoop, p.Epsi, g, mds);
							if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
								sigma = sign((xLoop[j])[dim-2]);
								ErgVor_nr(xLoop, j, mu, k, p.Epsi, sigma*g, erg);
							}
						}
						else if (poto==PotentialOptions::link)
							mdVlr_nr(j, mu, k, xLoop, p.Epsi, g, mds);
						else if (poto==PotentialOptions::exponential)
							mdVer_nr(j, mu, k, xLoop, p.Epsi, g, mds);
						else if (poto==PotentialOptions::dimreg)
							mdVdr_nr(j, mu, k, xLoop, p.Epsi, g, mds);
						else if (poto==PotentialOptions::thermal) {
							if (onlyselfreg) {
								mdVth2r_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
								if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
									sigma = sign((xLoop[j])[dim-2]);
									ErgVth2r_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*g, erg);
								}
							}
							else {
								mdVthr_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
								if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
									sigma = sign((xLoop[j])[dim-2]);
									ErgVthr_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*g, erg);
								}
							}
						}
						else if (poto==PotentialOptions::thermalDisjoint) {
							mdVthrDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
							if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
								sigma = sign((xLoop[j])[dim-2]);
								ErgVthrDisjoint_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*g, erg);
							}
						}
						else if (poto==PotentialOptions::thermalDisjointLR) {
							mdVthrDisjointLR_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
							if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
								sigma = sign((xLoop[j])[dim-2]);
								ErgVthrDisjointLR_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*g, erg);
							}
						}
						
						if (gaussian) {	
							if (!disjoint && poto!=PotentialOptions::thermal) {
								mdGaussian_nr(j, mu, k, xLoop, p.Epsi, repulsion_scale, mds);
								if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
									sigma = sign((xLoop[j])[dim-2]);
									ErgGaussian_nr(xLoop, j, mu, k, p.Epsi, sigma*repulsion_scale, erg);
								}
							}
							else if (disjoint && poto!=PotentialOptions::thermalDisjoint && poto!=PotentialOptions::thermalDisjointLR) {
								if (gaussianLR) {
									mdGaussianLRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, repulsion_scale, mds);
									if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
										sigma = sign((xLoop[j])[dim-2]);
										ErgGaussianLRDisjoint_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*repulsion_scale, erg);
									}
								}
								else {
									mdGaussianDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, repulsion_scale, mds);
									if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
										sigma = sign((xLoop[j])[dim-2]);
										ErgGaussianDisjoint_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*repulsion_scale, erg);
									}
								}
							}
							else if (poto==PotentialOptions::thermal) {
								//mdGaussian_nr(j, mu, k, xLoop, p.Epsi, repulsion_scale, mds);
								if (thermal2)
									mdGaussianThermal2_nr(j, mu, k, xLoop, beta, p.Epsi, repulsion_scale, mds);
								else 
									mdGaussian_nr(j, mu, k, xLoop, p.Epsi, repulsion_scale, mds);
								if (mu==(dim-1) && atCoord(xLoop,mu,0.0,j)) {
									sigma = sign((xLoop[j])[dim-2]);
									if (thermal2)
										ErgGaussianThermal2_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*repulsion_scale, erg);
									else
										ErgGaussian_nr(xLoop, j, mu, k, p.Epsi, sigma*repulsion_scale, erg);
								}
							}
							else if (poto==PotentialOptions::thermalDisjoint || poto==PotentialOptions::thermalDisjointLR) {
								if (gaussianLR) {
									if (thermal2)
										mdGaussianThermal2LRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, repulsion_scale, mds);
									else
										mdGaussianThermalLRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, repulsion_scale, mds);
									if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
										sigma = sign((xLoop[j])[dim-2]);
										if (thermal2)
											ErgGaussianThermal2LRDisjoint_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*repulsion_scale, erg);
										else
											ErgGaussianThermalLRDisjoint_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*repulsion_scale, erg);
									}
								}
								else {
									if (thermal2)
										mdGaussianThermal2Disjoint_nr(j, mu, k, xLoop, beta, p.Epsi, repulsion_scale, mds);
									else
										mdGaussianThermalDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, repulsion_scale, mds);
									if (mu==(dim-1) && atCoordDisjoint(xLoop,mu,0.0,j)) {
										sigma = sign((xLoop[j])[dim-2]);
										if (thermal2)
											ErgGaussianThermalDisjoint_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*repulsion_scale, erg);
										else
											ErgGaussianThermalDisjoint_nr(xLoop, j, mu, k, beta, p.Epsi, sigma*repulsion_scale, erg);
									}
								}
							}
						}
					}
				
					for (nu=0; nu<dim; nu++) {
					
						// free particle
						if (!disjoint) {
							if (kino==KineticOptions::saddle)
								ddsqrtS0_nr(j,mu,k,nu,xLoop,sqrt4s0,1.0,dds);
							else if (kino==KineticOptions::s0)
								ddS0_nr(j,mu,k,nu,xLoop,s0_scale,dds);
							else if (kino==KineticOptions::len)
								ddL_nr(j,mu,k,nu,xLoop,1.0,dds);
						}
						else {
							if (kino==KineticOptions::saddle)
								ddsqrtS0Disjoint_nr(j,mu,k,nu,xLoop,sqrt4s0,beta,1.0,dds);
							else if (kino==KineticOptions::s0)
								ddS0Disjoint_nr(j,mu,k,nu,xLoop,s0_scale,dds);
							else if (kino==KineticOptions::len)
								ddLDisjoint_nr(j,mu,k,nu,xLoop,beta,1.0,dds);
						}
						
						// external field
						if (!disjoint) {
							ddI0_nr(j,mu,k,nu,xLoop,mgb,dds);
							if (poto==PotentialOptions::external)
								ddIn_nr(j,mu,k,nu,xLoop,n,g,dds);
						}
						else {
							ddI0Disjoint_nr(j,mu,k,nu,xLoop,beta,mgb,dds);
							if (poto==PotentialOptions::externalDisjoint)
								ddInDisjoint_nr(j,mu,k,nu,xLoop,n,beta,g,dds);
						}
						
						
						if (!weak) {
							// dynamical field	
							if (poto==PotentialOptions::original) {
								ddVor_nr(j, mu, k, nu, xLoop, p.Epsi, g, dds);
							}
							else if (poto==PotentialOptions::link) {
								ddVlr_nr(j, mu, k, nu, xLoop, p.Epsi, g, dds);
							}
							else if (poto==PotentialOptions::exponential) {
								ddVer_nr(j, mu, k, nu, xLoop, p.Epsi, g, dds);
							}
							else if (poto==PotentialOptions::dimreg) {
								ddVdr_nr(j, mu, k, nu, xLoop, p.Epsi, g, dds);	
							}
							else if (poto==PotentialOptions::thermal) {
								if (onlyselfreg)
									ddVth2r_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
								else
									ddVthr_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
							}
							else if (poto==PotentialOptions::thermalDisjoint) {
								ddVthrDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
							}
							else if (poto==PotentialOptions::thermalDisjointLR) {
								ddVthrDisjointLR_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
							}
							else if (poto==PotentialOptions::nonrelDisjoint) {
								ddVnonrelrDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);		
							}
								
							// self-energy regularisation
							if (!disjoint) {
								if (gaussian && poto!=PotentialOptions::thermal) {
									ddGaussian_nr(j, mu, k, nu, xLoop, p.Epsi, repulsion_scale, dds);
								}
								else if (gaussian && poto==PotentialOptions::thermal) {
									if (thermal2)
										ddGaussianThermal2_nr(j, mu, k, nu, xLoop, beta, p.Epsi, repulsion_scale, dds);
									else
										ddGaussian_nr(j, mu, k, nu, xLoop, p.Epsi, repulsion_scale, dds);
									//ddGaussian_nr(j, mu, k, nu, xLoop, p.Epsi, repulsion_scale, dds);
								}
								else if (poto!=PotentialOptions::dimreg)
									ddL_nr(j,mu,k,nu,xLoop,dm,dds);
								else
									ddDistPow_nr(j,mu,k,nu,xLoop,p.Epsi,dim_reg_scale,dds);
							}
							else {
								if (gaussian && poto!=PotentialOptions::thermalDisjoint && poto!=PotentialOptions::thermalDisjointLR) {
									if (gaussianLR) {
										ddGaussianLRDisjoint_nr(j, mu, k, nu, xLoop, beta , p.Epsi, repulsion_scale, dds);
									}
									else {
										ddGaussianDisjoint_nr(j, mu, k, nu, xLoop, beta , p.Epsi, repulsion_scale, dds);
									}
								}
								else if (gaussian && (poto==PotentialOptions::thermalDisjoint || poto==PotentialOptions::thermalDisjointLR)) {
									if (gaussianLR) {
										if (thermal2)
											ddGaussianThermal2LRDisjoint_nr(j, mu, k, nu, xLoop, beta , p.Epsi, repulsion_scale, dds);
										else
											ddGaussianThermalLRDisjoint_nr(j, mu, k, nu, xLoop, beta , p.Epsi, repulsion_scale, dds);
									}
									else {
										if (thermal2)
											ddGaussianThermal2Disjoint_nr(j, mu, k, nu, xLoop, beta , p.Epsi, repulsion_scale, dds);
										else
											ddGaussianThermalDisjoint_nr(j, mu, k, nu, xLoop, beta , p.Epsi, repulsion_scale, dds);
									}
								}
								else
									ddLDisjoint_nr(j,mu,k,nu,xLoop,beta,dm,dds);
							}
						}								
					}
				}
			}		
		}
		
		// lagrange multiplier terms
		for (mu=0; mu<zm; mu++) {
			for (j=0; j<N; j++) {
				if ( (mu<(dim-1) && !fixall) || (mu==(dim-1) && !fixtlr && !fixall) ) {
					uint locj = j*dim+mu, locz = N*dim+mu;
					mds(locz) -= x[locj];
					mds(locj) -= x[locz];
				
					dds(locj,locz) += 1.0;
					dds(locz,locj) += 1.0;
				}
				else if (mu<=dim && fixtlr) {
					if (j<N/2 && mu==(dim-1)) {	// fixing average (dim-1) coord of LHS
						uint locj = j*dim+(dim-1), locz = N*dim+mu;
						mds(locz) -= x[locj];
						mds(locj) -= x[locz];
				
						dds(locj,locz) += 1.0;
						dds(locz,locj) += 1.0;
					}
					else if (j>=N/2&& mu==dim) {// fixing average (dim-1) coord of RHS
						uint locj = j*dim+(dim-1), locz = N*dim+mu;
						mds(locz) -= x[locj];
						mds(locj) -= x[locz];
				
						dds(locj,locz) += 1.0;
						dds(locz,locj) += 1.0;
					}
				}
				else if (mu<2*dim && fixall) {
					if (j<N/2 && !(mu%2)) {	// fixing average (uint)(mu/2) coord of LHS
						uint locj = j*dim+(uint)(mu/2), locz = N*dim+mu;
						mds(locz) -= x[locj];
						mds(locj) -= x[locz];
				
						dds(locj,locz) += 1.0;
						dds(locz,locj) += 1.0;
					}
					else if (j>=N/2 && mu%2) {// fixing average (uint)(mu/2) coord of RHS
						uint locj = j*dim+(uint)(mu/2), locz = N*dim+mu;
						mds(locz) -= x[locj];
						mds(locj) -= x[locz];
				
						dds(locj,locz) += 1.0;
						dds(locz,locj) += 1.0;
					}
				}
				else if (fixdz || fixodt) {
					uint mucf = (1+(uint)fixall)*dim+(uint)fixtlr;
					if ( (mu==mucf && j==(N/2-1) && fixodt) || (mu==(mucf+1) && j==(N-1) && fixodt && fixdislr) ) {
						// fixing relative heights of top and bottom points
						nu = dim-1;
						uint oj = oppNeigh(j,N);
						uint locj = j*dim+nu, locoj = oj*dim+nu, locz = N*dim+mu;
						mds(locz)  		-= x[locoj]-x[locj];
						mds(locoj)  	-= x[locz];
						mds(locj)  		-= -x[locz];								

						dds(locoj,locz)  	+= 1.0;
						dds(locz,locoj)  	+= 1.0;
						dds(locj,locz) 		+= -1.0;
						dds(locz,locj) 		+= -1.0;
						
					}
					if ( ( mu==mucf && j==(N/2-1) && fixdz) || ( mu==(mucf+1) && j==(N-1) && fixdz && fixdislr) ){
						uint nu = dim-2; // fixing dz=0 at top and bottom of RHS
						uint pj = (disjoint? posNeighDisjoint(j,N): posNeigh(j,N));
						uint locj = j*dim+nu, locpj = pj*dim+nu, locz = N*dim+mu;
						number ds = 1.0;///(number)N; // using 1/N makes mds large here
						mds(locz)  		-= (x[locpj]-x[locj])/ds;
						mds(locpj)  	-= x[locz]/ds;
						mds(locj)  		-= -x[locz]/ds;							

						dds(locpj,locz)  	+= 1.0/ds;
						dds(locz,locpj)  	+= 1.0/ds;
						dds(locj,locz) 		+= -1.0/ds;
						dds(locz,locj) 		+= -1.0/ds;
					}
				}
			}
		}
		
		if (!(P^=P0)) {
			// external momenta
			mdPX_nr(xLoop,N/2-1,P,-1.0,mds);
			mdPX_nr(xLoop,0,P,1.0,mds);
		
			// dynamical field cusp regularisation
			uint range = 2*p.Ng+1;
			uint left = N/2-1 - (int)p.Ng - 1;
			uint right = 0 + (int)(p.Ng>0)*(N-p.Ng);
			
			PseudoAngle(left, xLoop, 0.5, angle_neigh);
			left = (left==(N-1)? 0: left+1);
			
			for (k=0; k<range; k++) {			
				Angle(left, xLoop, 0.5/(number)range, gamma);
				Angle(right, xLoop, 0.5/(number)range, gamma);
				FGamma(left, xLoop, cusp_scale, fgamma);
				FGamma(right, xLoop, cusp_scale, fgamma);
				CuspCurvatureMax(left,xLoop,cc_scale,cc_max);
				CuspCurvatureMax(right,xLoop,cc_scale,cc_max);
				if (!weak) {
					mdFGamma_nr(xLoop,left,cusp_scale,mds);
					mdFGamma_nr(xLoop,right,cusp_scale,mds);
					ddFGamma_nr(xLoop,left,cusp_scale,dds);
					ddFGamma_nr(xLoop,right,cusp_scale,dds);
				}
				left = (left==(N-1)? 0: left+1);
				right = (right==(N-1)? 0: right+1);
			}
			
			PseudoAngle(right, xLoop, 0.5, angle_neigh);
		}
		
		// assigning scalar quantities
		vr = v;
		if (abs(p.Epsi)>MIN_NUMBER) {
			if (gaussian)
				vr += repulsion;
			else if (poto!=PotentialOptions::dimreg)
				vr += dm*len;
			else 
				vr += d_dim_reg;

			vr += (!(P^=P0)? cusp_scale*fgamma : 0.0);
		}
		if (kino==KineticOptions::saddle)
			kinetic = sqrt4s0;
		else if (kino==KineticOptions::s0)
			kinetic = s0_scale*s0;
		else if (kino==KineticOptions::len)
			kinetic = len;
		s = kinetic + i0;
		if (!weak) s += vr;
		if (!(P^=P0)) s -= Dot(xLoop[N/2-1]-xLoop[0],P);

		if (poto==PotentialOptions::thermal || disjoint)
			ergThermal = p.T*(kinetic+2.0*i0);
		else
			ergThermal = E;
		
/*----------------------------------------------------------------------------------------------------------------------------
	6 - some checks
----------------------------------------------------------------------------------------------------------------------------*/	

		// smoothness
		sm = Sm(xLoop,0,N/2-1);
		checkSm.add(sm);
		
		// check dx<<a, using the average
		number dx = len/(number)N;
		checkDX.add(dx/p.Epsi);
		
		// check ic<<1, a*kg<<1, dx*kg<<1, cc<<1, angle_neigh<<1
		checkICMax.add(ic_max);
		checkCCMax.add(cc_max);
		checkKgAMax.add(p.Epsi*kg_max);
		checkKgDxMax.add(dx*kg_max);
		checkAccAMax.add(p.Epsi*acc_max);
		checkStraight.add(angle_neigh);
		
		// conservation, Js
		number Js_mean = Js.sum()/(number)N;
		number Js_norm = Js.norm();
		Js = Js - Eigen::VectorXd::Constant(N,Js_mean);
		checkJs.add(Js.norm()/Js_norm);
		Js += Eigen::VectorXd::Constant(N,Js_mean);
					
		// check energies agree
		checkEAgree.add(2.0*abs(erg-ergThermal)/(abs(erg)+abs(ergThermal)));
				
		if (alltests) {
			// checking if dds is symmetric
			mat dds_asym(NT,NT);
			dds_asym = dds.transpose();
			dds_asym -= dds;
			dds_asym *= 0.5;
			number asym = dds_asym.norm();
			uint asymmaxx, asymmaxy, asymminx, asymminy;
			number asymmax = dds_asym.maxCoeff(&asymmaxx,&asymmaxy);
			number asymmin = dds_asym.minCoeff(&asymminx,&asymminy);
			if (-asymmin>asymmax) {
				asymmax = -asymmin;
				asymmaxx = asymminx;
				asymmaxy = asymminy;
			}
			checkSym.add(asym);
			checkSymMax.add(asymmax);
			checkSym.checkMessage();
			checkSymMax.checkMessage();
			cout << "position of max symmetric: (" << asymmaxx << "," << asymmaxy << ")/" << NT-1;
			cout << ", (j,k) = (" << (uint)(asymmaxx/dim) << "," << (uint)(asymmaxy/dim) << ")";
			cout << ", (mu,nu) = (" << asymmaxx%dim << "," << asymmaxy%dim << ")" << endl;
			
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
				cout << "negative eigenvalue = " << negEigenvalue << endl;
				cout << "analytic result     = " << analyticNegEigenvalue << endl;
				cout << "eigenvector test    = " << negEigenvectorTest << endl;
			}
			
			// checking if angle gamma agrees with weak coupling result
			number gamma_ratio = 0.0;
			if (abs(E)>MIN_NUMBER && abs(E)<2.0) {
				number gamma_free = 2.0*asin(E/2.0);
				gamma_ratio = gamma/gamma_free;
				checkGamma.add(gamma_ratio-1.0);
				checkGamma.checkMessage();
			}

			// check rotation and check mirror
			number xRotationTest = 0.0, xMirrorTest = 0.0;
			number mdsRotationTest = 0.0, mdsMirrorTest = 0.0;
			number offset = (runsCount==1? x[N*dim]: 0.00);
			number offsetT = (runsCount==1? x[N*dim]: 0.00); // used to be different if fixodt
			for (uint j=0; j<N/2; j++) {
				for (uint k=0; k<dim; k++) {
					if (k==(dim-2)) {
						xRotationTest += pow(x(dim*j+k)+x(dim*(N/2+j)+k),2);
						xMirrorTest += pow(x(dim*j+k)+x(dim*(N-1-j)+k),2);
						mdsRotationTest += pow(mds(dim*j+k)+mds(dim*(N/2+j)+k)+2.0*offset,2);
						mdsMirrorTest += pow(mds(dim*j+k)+mds(dim*(N-1-j)+k)+2.0*offset,2);
					}
					else if (k==(dim-1)) {
						if (!disjoint) {
							xRotationTest += pow(x(dim*j+k)+x(dim*(N/2+j)+k),2);
							xMirrorTest += pow(x(dim*j+k)-x(dim*(N-1-j)+k),2);
							mdsRotationTest += pow(mds(dim*j+k)+mds(dim*(N/2+j)+k)+2.0*offsetT,2);
							mdsMirrorTest += pow(mds(dim*j+k)-mds(dim*(N-1-j)+k),2);
						}
						else {
							xRotationTest += pow(mod<number>(x(dim*j+k)+x(dim*(N/2+j)+k),-beta/2.0,beta/2.0),2);
							xMirrorTest += pow(mod<number>(x(dim*j+k)-x(dim*(N-1-j)+k),-beta/2.0,beta/2.0),2);
							mdsRotationTest += pow(mod<number>(mds(dim*j+k)+mds(dim*(N/2+j)+k)+2.0*offsetT,-beta/2.0,beta/2.0),2);
							mdsMirrorTest += pow(mod<number>(mds(dim*j+k)-mds(dim*(N-1-j)+k),-beta/2.0,beta/2.0),2);
						}
					}
					else {
						xRotationTest += pow(x(dim*j+k)-x(dim*(N/2+j)+k),2);
						xMirrorTest += pow(x(dim*j+k)-x(dim*(N-1-j)+k),2);
						mdsRotationTest += pow(mds(dim*j+k)-mds(dim*(N/2+j)+k),2);
						mdsMirrorTest += pow(mds(dim*j+k)-mds(dim*(N-1-j)+k),2);
					}
				}
			}
			xRotationTest /= (xnorm/2.0);
			xMirrorTest /= (xnorm/2.0);
			mdsRotationTest /= (mds.squaredNorm()/2.0);
			mdsMirrorTest /= (mds.squaredNorm()/2.0);
			checkXRotation.add(sqrt(xRotationTest));
			checkXMirror.add(sqrt(xMirrorTest));
			checkMDSRotation.add(sqrt(mdsRotationTest));
			checkMDSMirror.add(sqrt(mdsMirrorTest));
		}		
	
/*----------------------------------------------------------------------------------------------------------------------------
	7 - print early 1
----------------------------------------------------------------------------------------------------------------------------*/	
		
		if (po!=PrintOptions::none) {
			Filename early = baseFolder+"data/temp/"+timenumber+"xEarly1_K_"+nts(p.K)+"_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(E)\
						+"_a_"+nts(p.Epsi)+".dat";
			if (abs(E)>MIN_NUMBER)
				(early.Extras).push_back(StringPair("mu",nts(p.Mu)));
			if (poto==PotentialOptions::thermal || disjoint)
				(early.Extras).push_back(StringPair("T",nts(p.T)));
			if (weak)
				(early.Extras).push_back(StringPair("weak","1"));
			if (poto!=PotentialOptions::original || gaussian)
				(early.Extras).push_back(potExtras);
			if (kino!=KineticOptions::saddle)
				(early.Extras).push_back(kinExtras);
				
			(early.Extras).push_back(StringPair("pl",nts(pl)));
			(early.Extras).push_back(StringPair("run",nts(runsCount)));
			
			if (po==PrintOptions::x || po==PrintOptions::all) {
				printAsLoop(early,dim,x,N*dim);
				//saveVectorAscii(early,x);
				printf("%12s%50s\n","x:",((string)early).c_str());
			}
			if (po==PrintOptions::mds || po==PrintOptions::all) {
				early.ID = "mdsEarly1";
				printAsLoop(early,dim,mds,N*dim);
				//saveVectorAscii(early,mds);
				printf("%12s%50s\n","mds:",((string)early).c_str());
			}
			if (po==PrintOptions::dds || po==PrintOptions::all) {
				early.ID = "ddsEarly1";
				//saveMatrixAscii(early,dds); // have stopped this because it just takes up too much space
				printf("%12s%50s\n","dds:",((string)early).c_str());
			}
			if (po==PrintOptions::curvature || po==PrintOptions::all) {
				early.ID = "curvatureEarly1";
				printAsLoop(early,dim,x,N*dim);
				saveVectorAsciiAppend(early,sc_vec);
				saveVectorAsciiAppend(early,kg_vec);
				saveVectorAsciiAppend(early,acc_vec);
				printf("%12s%50s\n","curvature:",((string)early).c_str());
			}
		}
		
/*----------------------------------------------------------------------------------------------------------------------------
	8 - solve for delta
----------------------------------------------------------------------------------------------------------------------------*/	
		
		// initializing delta
		number normx = x.norm();
		
		// solving for delta
		if (!pass) {				
			// solving for delta = DDS^{-1}*mdS
			delta = dds.partialPivLu().solve(mds);
			
		}
		else {
			delta = Eigen::VectorXd::Zero(NT);
		}
		
/*----------------------------------------------------------------------------------------------------------------------------
	9 - printing early 2 (delta), checking delta
----------------------------------------------------------------------------------------------------------------------------*/	

		// printing delta early
		if (po!=PrintOptions::none) {
			Filename early = baseFolder+"data/temp/"+timenumber+"deltaEarly2_K_"+nts(p.K)+"_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(E)\
							+"_a_"+nts(p.Epsi)+".dat";
			if (abs(E)>MIN_NUMBER)
				(early.Extras).push_back(StringPair("mu",nts(p.Mu)));
			if (poto==PotentialOptions::thermal || disjoint)
				(early.Extras).push_back(StringPair("T",nts(p.T)));
			if (weak)
				(early.Extras).push_back(StringPair("weak","1"));
			if (poto!=PotentialOptions::original || gaussian)
				(early.Extras).push_back(potExtras);
			if (kino!=KineticOptions::saddle)
				(early.Extras).push_back(kinExtras);
				
			(early.Extras).push_back(StringPair("pl",nts(pl)));
			(early.Extras).push_back(StringPair("run",nts(runsCount)));
			
			if (po==PrintOptions::delta || po==PrintOptions::all) {
					printAsLoop(early,dim,delta,N*dim);
				printf("%12s%50s\n","delta:",((string)early).c_str());
			}
		}		
		
		if (!pass && alltests) {
			// check rotation and check mirror of delta
			number deltaRotationTest = 0.0, deltaMirrorTest = 0.0;
			for (uint j=0; j<N/2; j++) {
				for (uint k=0; k<dim; k++) {
					if (k==(dim-2)) {
						deltaRotationTest += pow(delta(dim*j+k)+delta(dim*(N/2+j)+k),2);
						deltaMirrorTest += pow(delta(dim*j+k)+delta(dim*(N-1-j)+k),2);
					}
					else if (k==(dim-1)) {
						if (!disjoint) {
							deltaRotationTest += pow(delta(dim*j+k)+delta(dim*(N/2+j)+k),2);
							deltaMirrorTest += pow(delta(dim*j+k)-delta(dim*(N-1-j)+k),2);
						}
						else {
							deltaRotationTest += pow(mod<number>(delta(dim*j+k)+delta(dim*(N/2+j)+k),-beta/2.0,beta/2.0),2);
							deltaMirrorTest += pow(mod<number>(delta(dim*j+k)-delta(dim*(N-1-j)+k),-beta/2.0,beta/2.0),2);
						}
					}
					else {
						deltaRotationTest += pow(delta(dim*j+k)-delta(dim*(N/2+j)+k),2);
						deltaMirrorTest += pow(delta(dim*j+k)-delta(dim*(N-1-j)+k),2);
					}
				}
			}
			deltaRotationTest /= (delta.squaredNorm()/2.0);
			deltaMirrorTest /= (delta.squaredNorm()/2.0);
			checkDeltaRotation.add(sqrt(deltaRotationTest));
			checkDeltaMirror.add(sqrt(deltaMirrorTest));
			
			// for disjoint, checking force versus change
			if (disjoint) {
				number dsdz_RmL = 0.0;
				number deltaz_RmL = 0.0;
				uint mu = dim-2;
				for (uint j=0; j<N/2; j++) {
					dsdz_RmL 	+= -mds[dim*j+mu]-(-mds[dim*(j+N/2)+mu]);
					deltaz_RmL 	+= delta[dim*j+mu]-delta[dim*(j+N/2)+mu];
				}
				dsdz_RmL /= (number)N/2.0;
				deltaz_RmL /= (number)N/2.0;
				
				// printing summary
				string dzSummaryFile = "results/nr/nr_dz.csv";
				#define numDzSummary 19
				vector<string> dzSummary(numDzSummary);
				string dzSummary_array[numDzSummary] = {timenumber,\
											nts(pl),\
											nts(runsCount),\
											potExtras.second,\
											nts(p.K),\
											nts(pow(p.G,3)*p.B,16),\
											nts(p.Epsi,16),\
											nts(p.Mu,16),\
											nts(p.Lambda,16),\
											nts(E,16),\
											nts(p.T,16),\
											nts(s,16),\
											nts(zmax,16),\
											nts(dsdz_RmL,16),\
											nts(deltaz_RmL,16),\
											nts(checkSol.back(),16),\
											nts(checkDX.back(),16),\
											nts(checkDelta.back(),16),\
											nts(checkInv.back(),16)};								
				dzSummary.assign(dzSummary_array,dzSummary_array+numDzSummary);							
				saveVectorCsvAppend(dzSummaryFile,dzSummary);
				printf("%12s%50s\n","dz summary:",dzSummaryFile.c_str());
			}
		}

/*----------------------------------------------------------------------------------------------------------------------------
	10 - convergence checks
----------------------------------------------------------------------------------------------------------------------------*/	
	
		// calculating norms etc
		number normmds = mds.norm();
		uint maxmdspos, minmdspos;
		number maxdelta = 0.0;
		uint maxdeltapos = 0;
		number maxmds = mds.maxCoeff(&maxmdspos);
		number minmds = mds.minCoeff(&minmdspos);
		if (-minmds>maxmds) {
			maxmds = -minmds;
			maxmdspos = minmdspos;
		}
		number maxx = x.maxCoeff();
		number minx = x.minCoeff();
		if (-minx>maxx) {
			maxx = -minx;
		}
		number avgzm = (mds.tail(zm)).mean();
		number avgnzm = (mds.head(N)).mean();

		// adding to checks
		checkSol.add(normmds/normx);
		checkSolMax.add(maxmds/maxx);
		checkSolZM.add(avgzm/avgnzm);
		
		
		if (!pass) {
			number normdelta = delta.norm();
			uint mindeltapos;
			maxdelta = delta.maxCoeff(&maxdeltapos);
			number mindelta = delta.minCoeff(&mindeltapos);
			if (-mindelta>maxdelta) {
				maxdelta = -mindelta;
				maxdeltapos = mindeltapos;
			}
			checkDelta.add(normdelta/normx);
			
			// checking delta
			checkDelta.checkMessage();
			
			//independent check on whether calculation worked		
			number invError = (dds*delta - mds).norm();
			checkInv.add(invError);
			checkInv.checkMessage();
			
			if (!checkDelta.good() || !checkInv.good() || alltests) {
				printErrorInformation(x,"x",zm);
				cout << endl;
				printErrorInformation(mds,"mds",zm);
				cout << endl;
				printErrorInformation(delta,"delta",zm);
				cout << endl;
				printErrorInformation(dds,"dds");
				cout << endl;
				
				cout << "action info:" << endl;
				cout << "s:                      " << s               << endl;
				cout << "kinetic:                " << kinetic         << endl;
				cout << "i0:                     " << i0              << endl;
				cout << "vr:                     " << vr      << endl;
				if (!checkDelta.good() || !checkInv.good())
					passThrough = true;
			}
			
			//assigning values to x
			x += delta;
			
		}
	
		//printing tests to see convergence
		if (verbose) {
			if (runsCount==1) {
				printf("%4s%4s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","pl","run","len","i0","s","E","sol","solM","delta","a*acc_max","ic_max","cc_max","dx/a");
			}
			printf("%4i%4i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",pl,runsCount,len,i0,s,erg,\
			checkSol.back(),checkSolMax.back(),checkDelta.back(),checkAccAMax.back(),\
				checkICMax.back(),checkCCMax.back(),checkDX.back());
		}
		if (alltests) {
			checkSolZM.checkMessage();
			checkKgAMax.checkMessage();
			checkStraight.checkMessage();
			checkXRotation.checkMessage();
			checkXMirror.checkMessage();
			checkMDSRotation.checkMessage();
			checkMDSMirror.checkMessage();
			checkDeltaRotation.checkMessage();
			checkDeltaMirror.checkMessage();
			cout << "avg mds:               " << mds.mean() << endl;
			cout << "max mds:               " << maxmds << endl;
			cout << "position of max mds:   " << maxmdspos << "/" << NT-1 << ", j = " << (uint)(maxmdspos/dim);
			cout << ", mu = " << maxmdspos%dim << endl;
			cout << "avg delta:             " << delta.mean() << endl;
			if (!pass) {
				cout << "max delta:             " << maxdelta << endl;
				cout << "position of max delta: " << maxdeltapos << "/" << NT-1 << ", j = " << (uint)(maxdeltapos/dim);
				cout << ", mu = " << maxdeltapos%dim << endl;
			}
		}
		if (conservation || alltests || sometests) {
			checkJs.checkMessage();
			checkEAgree.checkMessage();
			if (!checkEAgree.good()) {
				cout << "erg = " << erg << ", ergThermal = " << ergThermal;
				if (p.T<0.5) {
					cout << ", ergLowTemp = " << 2.0*pow(PI,4)*(pow(p.G,3)*p.B)*pow(p.T,5)/45.0 << endl;
				}
				else {
					cout << ", ergHighTemp = " << 2.0*(1.0-sqrt((pow(p.G,3)*p.B)/4.0/PI)) << endl;
				}
			}
			string consFile = baseFolder+"data/temp/"+timenumber+"Js_pl_"+nts(pl)+"_run_"+nts(runsCount)+".dat";
			saveVectorAscii(consFile,Js);
			printf("%12s%50s\n","Js       :",consFile.c_str());
		}
	
	}
	
/*----------------------------------------------------------------------------------------------------------------------------
	11 - printing results
----------------------------------------------------------------------------------------------------------------------------*/	

	// eigenvalues, if required
	if (eigen) {
		mat dds_wlm;
		if (includeLagrangeMultipliers)
			dds_wlm = dds;
		else
			dds_wlm = dds.block(0,0,dim*N,dim*N); // dds without Lagrange multipliers
		number eigenTol = MIN_NUMBER*pow(dim*N,2);
		number cos = 0.0;
		uint negEigs = 0;
		uint zeroEigs = 0;
		uint numEigs = 8*dim;
		cout << "calculating eigendecomposition of dds..." << endl;
		Eigen::SelfAdjointEigenSolver<mat> eigensolver(dds_wlm);
		if (eigensolver.info()!=Eigen::Success) abort();
		Filename eigenFile = baseFolder+"data/nr/eigenvalues/dim_"+nts(dim)+"/K_"+nts(p.K)+"/"+timenumber+"eigenvalues_pl_"+nts(pl)\
				+"_run_"+nts(runsCount)+".dat";
		saveVectorBinary(eigenFile,eigensolver.eigenvalues());
		printf("%12s%50s\n","eigenvalues:",((string)eigenFile).c_str());
		cout << "first " << numEigs << " eigenvalues and their dot products with mds are: " << endl;
		/////////////////////////////////////////
		/*cout << "quick eigs test:" << endl;
		cout << "<v0|v0> = " << ((eigensolver.eigenvectors()).col(0)).norm() << endl;
		cout << "<v0|DDS|v0> = " << ((eigensolver.eigenvectors()).col(0)).dot(dds_wlm*((eigensolver.eigenvectors()).col(0))) << endl;*/
		////////////////////////////////////////
		for (uint j=0; j<numEigs; j++) {
			if ((eigensolver.eigenvalues())[j]<-eigenTol)
				negEigs++;
			if (abs((eigensolver.eigenvalues())[j])<eigenTol)
				zeroEigs++;
			if (includeLagrangeMultipliers)
				cos = mds.dot((Eigen::VectorXd)(eigensolver.eigenvectors()).col(j))/mds.norm();
			else
				cos = (mds.head(dim*N)).dot((Eigen::VectorXd)(eigensolver.eigenvectors()).col(j))/(mds.head(dim*N)).norm();
			cout << (eigensolver.eigenvalues())[j] << " " << cos << endl;
			eigenFile.ID = "eigenvector"+nts(j);
			saveVectorBinary(eigenFile,(Eigen::VectorXd)((eigensolver.eigenvectors()).col(j)));
		}
		cout << negEigs << " negative eigenvalues found, less than " << -eigenTol << endl;
		cout << zeroEigs << " zero eigenvalues found, with absolute value less than " << eigenTol << endl;
		printf("%12s%50s\n","eigenvectors:",((string)eigenFile).c_str());
		
		// printing eigenvalues to file
		// id
		vector<string> idEigs(idSizeResults);
		idEigs[0] = timenumber;
		idEigs[1] = nts(pl);
		idEigs[2] = potExtras.second;
		
		// actual results
		vector<number> datumEigs(11+numEigs);
		datumEigs[0] = s;
		datumEigs[1] = ergThermal;
		datumEigs[2] = zmax;
		datumEigs[3] = zmin;
		datumEigs[4] = tmax;
		datumEigs[5] = checkSol.back();
		datumEigs[6] = checkDX.back();
		datumEigs[7] = checkICMax.back();
		datumEigs[8] = checkAccAMax.back();
		datumEigs[9] = checkCCMax.back();
		datumEigs[10] = checkStraight.back();
		
		for (uint j=0; j<numEigs; j++)
			datumEigs[11+j] = (eigensolver.eigenvalues())[j];
		
		// saving
		string eigsFile = "results/nr/nr_eigs8.csv";
		NewtonRaphsonDatum eigs(idEigs,p,datumEigs);
		eigs.save(eigsFile);
		if (verbose)
			printf("%12s%24s\n","eigenvalues:",eigsFile.c_str());
	}
	
	// curvature, if required
	if (curvature || alltests) {
		Filename file = baseFolder+"data/temp/"+timenumber+"xCurvature_K_"+nts(p.K)+"_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(E)\
						+"_a_"+nts(p.Epsi)+".dat";
		if (abs(E)>MIN_NUMBER)
			(file.Extras).push_back(StringPair("mu",nts(p.Mu)));
		if (poto==PotentialOptions::thermal || disjoint)
			(file.Extras).push_back(StringPair("T",nts(p.T)));
		printAsLoop(file,dim,x,N*dim);
		saveVectorAsciiAppend(file,sc_vec);
		saveVectorAsciiAppend(file,kg_vec);
		saveVectorAsciiAppend(file,acc_vec);
		printf("%12s%50s\n","curvature:",((string)file).c_str());
	}

	//stopping clock
	time = clock() - time;
	number realtime = time/1000000.0;
	
	// results to compare to
	number s_cf = PI/p.G/p.B;
	if (!weak) s_cf -= p.G*p.G/4.0;
	else if (abs(E)>MIN_NUMBER && abs(E)<2.0) s_cf -= (2.0/p.G/p.B)*( asin(E/2.0) + (E/2.0)*sqrt(1.0-pow(E/2.0,2)) );
		
	// printing results to terminal
	printf("\n");
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%12s%12s%14s%14s%14s%14s%14s\n","pl","runs","time","K","G","B","Ng","a","mu","E","T","len",\
		"vr","s","acca","dxa");
	printf("%8i%8i%8.3g%8i%8.4g%8.4g%8i%8.4g%8.4g%12.5g%12.5g%14.5g%14.5g%14.5g%14.5g%14.5g\n",\
		pl,runsCount,realtime,p.K,p.G,p.B,p.Ng,p.Epsi,p.Mu,erg,p.T,len,vr,s,checkAccAMax.back(),checkDX.back());
	printf("\n");
	
	if (checkDelta.good() && checkSol.good() && checkSolMax.good()) {		
		// printing loop to file
		Filename loopRes;
		if (poto==PotentialOptions::thermal || disjoint)
			loopRes = filenameThermalNR<dim>(p,baseFolder);
		else
			loopRes = filenameLoopNR<dim>(p,baseFolder);
		if (weak)
			(loopRes.Extras).push_back(StringPair("weak","1"));
		if (poto!=PotentialOptions::original || gaussian)
			(loopRes.Extras).push_back(potExtras);
		if (kino!=KineticOptions::saddle)
			(loopRes.Extras).push_back(kinExtras);
		saveVectorBinary(loopRes,x);
		if (verbose)
			printf("%12s%50s\n","x:",((string)loopRes).c_str());
	}
	
	if ((checkDelta.good() && checkSol.good() && checkSolMax.good()) || pass) {
		// printing good results to file
		
		// id
		vector<string> idResult(idSizeResults);
		idResult[0] = timenumber;
		idResult[1] = nts(pl);
		idResult[2] = potExtras.second;
		
		// actual results
		vector<number> datumResult(datumSizeResults);
		datumResult[0] = s;
		datumResult[1] = erg;
		datumResult[2] = ergThermal;
		datumResult[3] = gamma;
		datumResult[4] = len;
		datumResult[5] = kinetic;
		datumResult[6] = i0;
		datumResult[7] = vr;
		datumResult[8] = zmax;
		datumResult[9] = zmin;
		datumResult[10] = tmax;
		datumResult[11] = checkSol.back();
		datumResult[12] = checkDX.back();
		datumResult[13] = checkICMax.back();
		datumResult[14] = checkAccAMax.back();
		datumResult[15] = checkCCMax.back();
		datumResult[16] = checkStraight.back();
		
		// saving
		NewtonRaphsonDatum result(idResult,p,datumResult);
		result.save(resultsFile);
		if (verbose)
			printf("%12s%24s\n","results:",resultsFile.c_str());
	}
	else {
		// printing error results to file	
		
		// id
		vector<string> idError(idSizeErrors);
		idError[0] = timenumber;
		idError[1] = nts(pl);
		idError[2] = nts(runsCount);
		idError[3] = potExtras.second;
		
		// actual results
		vector<number> datumError(datumSizeErrors);
		datumError[0] = s;
		datumError[1] = erg;
		datumError[2] = ergThermal;
		datumError[3] = len;
		datumError[4] = checkSol.back();
		datumError[5] = checkSolMax.back();
		datumError[6] = checkDelta.back();
		datumError[7] = checkDX.back();
		datumError[8] = checkICMax.back();
		datumError[9] = checkAccAMax.back();
		datumError[10] = checkCCMax.back();
		datumError[11] = checkStraight.back();
		datumError[12] = checkJs.back();
		
		// saving
		NewtonRaphsonDatum error(idError,p,datumError);
		error.save(errorsFile);
		if (verbose)
			printf("%12s%24s\n","results:",errorsFile.c_str());
	}
	
	// printing stepper results
	if (stepargv!=StepperArgv::none) {
		// adding result to stepper
		number F = 0.0;
		if(stepargv==StepperArgv::action)
			F = s;
		else if(stepargv==StepperArgv::entropy)
			F = s - ergThermal*p.T;
		else {
			cerr << "nrmain_fn error: stepargv, " << stepargv << ", not recognized" << endl;
			return 1;
		}
		stepper.addResult(F);
		
		// printing step
		stepper.save(stepperOutputsFile);
		if (verbose)
			printf("%12s%24s\n","stepper outputs:",((string)stepperOutputsFile).c_str());
		
		// stepping
		stepper.step();
	
	}

	// printing extras to ascii files
	if (po!=PrintOptions::none) {
		Filename file = baseFolder+"data/temp/"+timenumber+"xEnd_K_"+nts(p.K)+"_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(E)\
							+"_a_"+nts(p.Epsi)+".dat";
			if (abs(E)>MIN_NUMBER)
				(file.Extras).push_back(StringPair("mu",nts(p.Mu)));
			if (poto==PotentialOptions::thermal || disjoint)
				(file.Extras).push_back(StringPair("T",nts(p.T)));
			if (weak)
				(file.Extras).push_back(StringPair("weak","1"));
			if (poto!=PotentialOptions::original || gaussian)
				(file.Extras).push_back(potExtras);
			if (kino!=KineticOptions::saddle)
				(file.Extras).push_back(kinExtras);
				
		if (po==PrintOptions::x || po==PrintOptions::all) {
			printAsLoop(file,dim,x,N*dim);
			printf("%12s%50s\n","x:",((string)file).c_str());
		}
		else if (po==PrintOptions::mds || po==PrintOptions::all) {
			file.ID = "mdsEnd";
			printAsLoop(file,dim,x,N*dim);
			printf("%12s%50s\n","mds:",((string)file).c_str());
		}
	}
}

return 0;
}
