/*
	selfintreg
		quick program to load binary files containing a vector of doubles and print ascii ones
*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include "simple.h"

using namespace std;

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

uint min = 2, max = 10;
number a = 0.0;

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("min")==0) min = stn<uint>(argv[2*j+2]);
		else if (id.compare("max")==0) max = stn<uint>(argv[2*j+2]);
		else if (id.compare("a")==0) a = stn<number>(argv[2*j+2]);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. calculating and printing
----------------------------------------------------------------------------------------------------------------------------*/

number logFr, theta2, Fr2; //Fr = 0, theta
uint N, m;

//printf("%20s%20s%20s%20s\n","K","N","Fr","log(Fr)");
for (uint n=min; n<=max; n++) {
	N = pow(2,n);
	//Fr = 0.0;
	Fr2 = 0.0;
	for (uint i=2; i<=N; i++) {
		m = i-1;
		theta2 = pi*(number)m/(number)N;
		Fr2 += (number)(N-m)*gsl_sf_cos(theta2)/(1.0-gsl_sf_cos(theta2)+a*a);
		/*for (uint j=1; j<i; j++) {
			theta = pi*(number)(i-j)/(number)N;
			Fr += gsl_sf_cos(theta)/(1.0-gsl_sf_cos(theta)+a*a);
		}*/
	}
	//Fr *= 2.0*(1.0-gsl_sf_cos(pi/(number)N));
	Fr2 *= 2.0*(1.0-gsl_sf_cos(pi/(number)N));
	logFr = gsl_sf_log(Fr2)/gsl_sf_log(2.0);
	printf("%20i%20i%20g%20g\n",n,N,Fr2,logFr);
}

return 0;
}
