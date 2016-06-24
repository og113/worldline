/*
	quick test of parameters
*/

#include <iostream>
#include <fstream>
#include <string>
#include "simple.h"
#include "parameters.h"

using namespace std;

int main(int argc, char** argv) {

// data to print
string inputsFile = "inputs4";
string outputsFile = "temp/inputs";
bool verbose = true;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("inputs")==0 || id.compare("in")==0 || id.compare("i")==0) inputsFile = (string)argv[2*j+2];
		else if (id.compare("outputs")==0 || id.compare("out")==0 || id.compare("o")==0) outputsFile = (string)argv[2*j+2];
		else if (id.compare("verbose")==0) verbose = (stn<uint>(argv[2*j+2])!=0);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

if (verbose) {
	printf("%12s%24s\n","inputs:  ",inputsFile.c_str());
}

// loading parameters
ParametersRange pr, prOut;
pr.load(inputsFile);
Parameters p = pr.Min;

// parameter loops
//uint Npl = pr.totalSteps();
Parameters::Label labelT = Parameters::t;

uint Npl = 1; // number of parameter loops

if (pr.totalSteps()>1) {
	cout << "N.B. highTempParameters can only step over T" << endl;
}
if (pr.toStep(labelT)) {
	Npl = pr.totalSteps();
}

// starting loop
for (uint pl=0; pl<Npl; pl++) {

	bool pass = false;
	while (!pass) {

		// stepping parameters
		if (pr.toStep(labelT) && pl>0) {
			p = pr.position(pl);
		}
	
		string outputsFile_pl = outputsFile + "_T_"+nts(p.T);
	
		prOut.Min = p;
		prOut.Max = p;

		number gg = ((prOut.Min).Lambda>1? (prOut.Min).Lambda: 3.0);
		number T = (prOut.Min).T;
		number beta = 1.0/T;

		// steps
		uint steps = ((prOut.Min).Ng>1? (prOut.Min).Ng: 5);
		fill((prOut.Steps).begin(),(prOut.Steps).end(),0.0);
		(prOut.Steps)[7] = steps;
		(prOut.Steps)[9] = steps;

		// kappa
		number kappaMax = 1.0/(4.0*pow(PI,3)*pow(T,4)); // so that nonrelativistic solution exists
		number c = 2.22341;
		number B = c/(2.0*sqrt(2.0)*pow(PI,1.25));
		number deltaKappa = B*kappaMax/gg; // deltaKappa << B*kappaMax, so that dr<<r (harder than dr<<beta if sqrt(kappa/4.0/PI)<beta)
		number kappaMin = pow(beta,6)/pow(PI,3);
		number kappaMinSoft = kappaMax - deltaKappa;
		if (kappaMax<kappaMin) {
			cerr << "highTempParamters Error: on pl=" << pl << ", kappaMax(" << kappaMax << ") < kappaMin(" << kappaMin << ")" << endl;
			pass = true;
		}
		else if (kappaMin>kappaMinSoft)
			kappaMinSoft = kappaMin;
		(prOut.Min).B = kappaMinSoft;
		(prOut.Max).B = kappaMax;
		(prOut.Min).G = 1.0;
		(prOut.Max).G = 1.0;

		// a
		number N = pow(2.0,(prOut.Min).K);
		number dl = (2.0/T)/N;
		number aMin = gg*dl;
		number r0Max = sqrt(kappaMax/4.0/PI);
		number r0Min = sqrt(kappaMin/4.0/PI);
		number drMax = (pow(PI,1.0/8.0)/pow(2.0,5.0/4.0)/sqrt(c))*pow(kappaMin,3.0/8.0)*sqrt(deltaKappa);
		//number drMin = 0.0;
		number aMax = r0Max/gg; // so a<<r (a<<dr is harder but not so obviously necessary)
		if (aMax<aMin) {
			cerr << "highTempParamters Error: on pl=" << pl << ", aMax(" << aMax << ")<aMin(" << aMin << ")" << endl;
			pass = true;
		}
		(prOut.Min).Epsi = aMin;
		(prOut.Max).Epsi = aMax;

		// printing parameters to terminal
		if (verbose && !pass) {
			printf("\n");
			printf("%14s%14s%14s%14s%14s%14s%14s%14s%14s\n","kappaMin","kappaMax","aMin","aMax","beta","r0Min","r0Max","drMax","dl");
			printf("%14.5g%14.5g%14.5g%14.5g%14.5g%14.5g%14.5g%14.5g%14.5g\n",\
				kappaMin,kappaMax,aMin,aMax,1.0/T,r0Min,r0Max,drMax,dl);
			printf("\n");
		}

		if (!pass)
			prOut.save(outputsFile_pl);
		if (verbose && !pass)
			printf("%12s%24s\n","outputs: ",outputsFile_pl.c_str());

		pass = true;
		}
	}
}
