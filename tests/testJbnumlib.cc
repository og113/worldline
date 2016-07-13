// test of jbnumlib which i found along with the package called chiron (on arxiv)
#include <iostream>
#include <iomanip>
#include <vector>
#include <jbnumlib.h>
#include "simple.h"

using namespace std;

// Gaussian
number Gaussian(const number& r, const number& t, const number& a) {
	return exp(-(pow(r,2)+pow(t,2))/pow(a,2));
}


// GThermal
number GThermal(const number& r, const number& t, const number& beta, const number& a) {
	return (a*sqrt(PI)*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/(exp(pow(r,2)/pow(a,2))*beta);
}

int main() {

cout << "test of jbnumlib:" << endl;

cout << setw(20) << "u" << setw(20) << "q";
cout << setw(20) << "jbdtheta3(u,q)" << setw(20) << "jbderivutheta3(u,q)" << setw(20) << "jbderiv2utheta3(u,q)" << endl;

uint N = 4;
number u, q;
cout << setprecision(10);
for (uint j=0; j<N; j++) {
	u = 2.0*j/(number)N;
	for (uint k=0; k<N; k++) {
		q = k/(number)N;
	
		cout << setw(20) << u << setw(20) << q;
		cout << setw(20) << jbdtheta3(u,q) << setw(20) << jbderivutheta3(u,q) << setw(20) << jbderiv2utheta3(u,q) << endl;
	}
}

cout << setw(20) << "u" << setw(20) << "q";
cout << setw(20) << "jbdtheta3(u,q)" << setw(20) << "jbderivutheta3(u,q)" << setw(20) << "jbderiv2utheta3(u,q)" << endl;

cout << "comparison of Gaussian and GThermal:" << endl;

number beta = 1.0, a = 0.3;
cout << "beta = " << beta << endl;
cout << "a = " << a << endl;

cout << setw(20) << "r" << setw(20) << "t";
cout << setw(20) << "Gaussian(r,t,a)" << setw(20) << "GThermal(r,t,beta,a)" << setw(20) << "difference" << endl;

cout << setprecision(10);
number r, t;
for (uint j=0; j<N; j++) {
	r = 0.01*j/(number)N;
	for (uint k=0; k<N; k++) {
		t = 2.0*k/(number)N;
	
		cout << setw(20) << r << setw(20) << t;
		cout << setw(20) << Gaussian(r,t,a) << setw(20) << GThermal(r,t,beta,a) << setw(20) << Gaussian(r,t,a)-GThermal(r,t,beta,a) << endl;
	}
}

return 0;
}
