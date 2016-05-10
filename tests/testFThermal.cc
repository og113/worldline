#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include "simple.h"

// F
number F(const number& r, const number& t, const number& beta, const number& a) {
	return -1.0/(pow(2.0*PI,2)*(pow(a,2)+pow(r,2)+pow(t,2)));
}

// sinch
static number sinch(const number& x) {
	return (abs(x)<MIN_NUMBER? 1.0: sinh(x)/x);
}

// coshsinhc
static number coshsinhc(const number& x) {
	return (abs(x)<MIN_NUMBER? 0.0: cosh(x)/x - sinh(x)/pow(x,2));
}

// coshsinhcc
static number coshsinhcc(const number& x) {
	return (abs(x)<MIN_NUMBER? 1.0/3.0: cosh(x)/pow(x,2) - sinh(x)/pow(x,3));
}

// FThermal
number FThermal(const number& r, const number& t, const number& beta, const number& a) {
	return sinch((2.0*r*PI)/beta)/(2.0*pow(beta,2)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*r*PI)/beta)));
}

// DFThermalDr
number DFThermalDr(const number& r, const number& t, const number& beta, const number& a) {
	return (PI*coshsinhc((2.0*PI*r)/beta))/\
 (pow(beta,3)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta))) \
 + (2.0*pow(PI,2)*r*pow(sinch((2.0*PI*r)/beta),2))/\
 (pow(beta,4)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
}

// DFThermalDt
number DFThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	return (PI*sin((2.0*PI*t)/beta)*sinch((2.0*PI*r)/beta))/\
 (pow(beta,3)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
}

// DDFThermalDrDr
number DDFThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*pow(PI,2)*(-((2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta))*
 (3.0*cosh((2.0*PI*r)/beta) - 2.0*sinch((2.0*PI*r)/beta))*sinch((2.0*PI*r)/beta)) \
 + 8.0*pow(PI,2)*pow(r,2)*pow(sinch((2.0*PI*r)/beta),3) \
 + (pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),2)*
 (-2.0*coshsinhcc((2.0*PI*r)/beta) + sinch((2.0*PI*r)/beta)))/pow(beta,2)))/\
 pow(-2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cos((2.0*PI*t)/beta) - pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}

// DDFThermalDtDr
number DDFThermalDtDr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*pow(PI,2)*sin((2.0*PI*t)/beta)*((2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta))*
 coshsinhc((2.0*PI*r)/beta) - 4.0*PI*r*beta*pow(sinch((2.0*PI*r)/beta),2)))/\
 pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}

// DDFThermalDtDt
number DDFThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	return  (pow(PI,2)*(-3.0*pow(beta,2) + pow(beta,2)*cos((4.0*PI*t)/beta) \
 + 2.0*cos((2.0*PI*t)/beta)*(2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cosh((2.0*PI*r)/beta)))*sinch((2.0*PI*r)/beta))/\
 pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}

int main() {

cout << "test FThermal" << endl;

uint N=20;
number beta = 1.0e3;
number a = 0.1;
number r = 0.0, t = 0.0;
number rmax = 2.0, tmax = 2.0;

cout << setw(20) << "F" << setw(20) << "FThermal" << setw(20) << "DFThermalDr" << setw(20) << "DDFThermalDrDr" << endl;
for (uint j=0; j<N; j++) {
	r = rmax*j/(number)(N-1.0);
	t = tmax*j/(number)(N-1.0);
	cout << setw(20) << F(r,t,beta,a) << setw(20) << FThermal(r,t,beta,a) << setw(20) << DFThermalDr(r,t,beta,a);
	cout <<setw(20) << DDFThermalDrDr(r,t,beta,a) << endl;
}

return 0;
}
