/*
	testdF
		writing a function to take as an argument the scalar functions dFdt etc and then to spit our mdF as a vector
*/

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "analysis.h"
#include "check.h"
#include "folder.h"
#include "genloop.h"
#include "nrloop.h"
#include "print.h"
#include "simple.h"

using namespace std;

typedef number(*ThermalFunction)(const number&, const number&, const number&, const number&);

// FThermal
number FThermal(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return  sinh((2.0*PI*ra)/beta)/(4.0*PI*ra*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))); //
}

// DFThermalDrOnr
number DFThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return cosh((2.0*PI*ra)/beta)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 - sinh((2.0*PI*ra)/beta)/\
 (4.0*PI*pow(ra,3)*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + pow(sinh((2.0*PI*ra)/beta),2)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2));
 }

// DFThermalDt
number DFThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (sin((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (2.0*ra*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2));
}
 
 // mdVthr_nr
template<uint Dim>
void mdFth_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, \
					 	const number& beta, const number& a, const number& f, vec& v);

// mdFth_nr
template <> void mdFth_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
			ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, \
			const number& beta, const number& a, const number& f, vec& v) {
		
	number dimension = 4;	
	number res = 0.0;
	
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DX(l,j,i,dimension-1); // checked order
		number FThermal_ij = Fth(r_ij,t_ij,beta,a);
		number DFThermalDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
		number DFThermalDt_ij = dFthdt(r_ij,t_ij,beta,a);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		
		res += 2.0*FThermal_ij*(-DX(l,i,mu));
		if (mu<(dimension-1))
			res += 2.0*DFThermalDrOnr_ij*DX(l,j,i,mu)*T_ij;
		else
			res += 2.0*DFThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DX(l,mj,i,dimension-1); // checked order
		number FThermal_imj = Fth(r_imj,t_imj,beta,a);
		res +=  -2.0*FThermal_imj*(-DX(l,i,mu)); //
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (i==j)
		res += (-1.0/pow(2.0*PI,2))*2.0*(l[j])[mu]/a/a; 
	if (i==mj)
		res += (-1.0/pow(2.0*PI,2))*(-(l[mj])[mu]/a/a);
	if (i==pj)
		res += (-1.0/pow(2.0*PI,2))*(-(l[pj])[mu]/a/a);
		
	v[j*4+mu] += -f*(-pow(2.0*PI,2))*res;
}

int main() {

cout << "testdF:" << endl;

// inputs
string inputsFile = "inputs4";
cout << "loading parameters from " << inputsFile << endl;
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file: " << inputsFile << endl;
	return 1;
}

#define dim 4
uint N = pow(2,p.K);
uint zm = dim+1;
uint NT = N*dim+zm;

// defining vector and matrix quantities
vec x(N*dim);
vec mds(NT);

// defining xLoop
Loop<dim> xLoop(p.K,0);

string loadFile = "data/nr/loops/dim_4/K_11/loop_kappa_0.1_E_0_a_0.02_T_0.02_pot_8.dat";
cout << "loading x from " << loadFile << endl;

// loading x
loadVectorBinary(loadFile,x);
if (x.size()<NT) {
	x.conservativeResize(NT);
	for (uint mu=0; mu<zm; mu++)
		x[N*dim+mu] = 1.0e-4;
}
else if (x.size()>NT)
	x.conservativeResize(NT);
	
// loading x to xLoop - messier than it should be (should work with either a vec or a Loop really)
vectorToLoop(x,xLoop);

number beta = 1.0/p.T;
number g = p.G*p.G/8.0/PI;

clock_t time;
time = clock();

cout << "norm(mds) before: " << mds.norm() << endl;

// bulk
for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dim; mu++) {
		for (uint k=0; k<N; k++) {
			mdFth_nr(j, mu, k, xLoop, (ThermalFunction)&FThermal, (ThermalFunction)&DFThermalDrOnr, (ThermalFunction)&DFThermalDt, beta, p.Epsi, g, mds);
		}
	}
}

cout << "norm(mds) middle: " << mds.norm() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;

// bulk
for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dim; mu++) {
		for (uint k=0; k<N; k++) {
			mdVthr_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
		}
	}
}

cout << "norm(mds) after: " << mds.norm() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;

cout << "before and after should be zero, middle nonzero." << endl;

return 0;
}
