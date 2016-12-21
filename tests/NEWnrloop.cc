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
#include "nrloop2.h"
#include "print.h"
#include "simple.h"

using namespace std;

#define dim 4

int main() {

cout << "==========================================================================================================" << endl << endl;

cout << "test new nrloop:" << endl;

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

bool disjoint = false;
bool gaussian = true;

// some scalar quantities
uint N = pow(2,p.K);
uint zm = dim+1;
uint NT = N*dim+zm;
number beta = 1.0/p.T;
number g = p.G*p.G/8.0/PI;
number v = 0.0, repulsion = 0.0;

// defining vector and matrix quantities
vec x(N*dim);
vec mds(NT);
mat dds(NT,NT);

// defining xLoop
Loop<dim> xLoop(p.K,0);

// file to load x from
string loadFile = "data/nr/loops/dim_4/K_6/loop_kappa_1_E_0_a_0.1_T_0.02_pot_8.dat";
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
	
// loading x to xLoop
vectorToLoop(x,xLoop);

clock_t time;

cout << "-------------------------------------------------------------------------------------------------------------------------" << endl;
cout << "dim = "<< dim << endl;

cout << "v before = " << v << endl;
cout << "repulsion before = " << repulsion << endl;
cout << "norm(mds) before: " << mds.norm() << endl;
cout << "norm(dds) before: " << dds.norm() << endl;

time = clock();

// new ones
for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dim; mu++) {
		for (uint k=0; k<N; k++) {
			if (mu==0) {
				NEWVor(j,k,xLoop,p.Epsi,g,v);
				//NEWVer(j,k,xLoop,p.Epsi,g,v);
				//NEWVdr(j,k,xLoop,p.Epsi,g,v);
				NEWVthr(j,k,xLoop,beta,p.Epsi,g,v);
				NEWGaussian(j, k, xLoop, p.Epsi, g, repulsion);
				if (disjoint) {
					NEWVthrDisjoint(j,k,xLoop,beta,p.Epsi,g,v);
					NEWGaussianDisjoint(j, k, xLoop, beta, p.Epsi, g, repulsion);
					NEWGaussianLRDisjoint(j, k, xLoop, beta, p.Epsi, g, repulsion);
					NEWGaussianThermal(j, k, xLoop, beta, p.Epsi, g, repulsion);
					NEWGaussianThermal2(j, k, xLoop, beta, p.Epsi, g, repulsion);
					NEWGaussianThermalDisjoint(j, k, xLoop, beta, p.Epsi, g, repulsion);
					NEWGaussianThermal2Disjoint(j, k, xLoop, beta, p.Epsi, g, repulsion);
					NEWGaussianThermalLRDisjoint(j, k, xLoop, beta, p.Epsi, g, repulsion);
					NEWGaussianThermal2LRDisjoint(j, k, xLoop, beta, p.Epsi, g, repulsion);
				}
			}
			NEWmdVor_nr(j, mu, k, xLoop, p.Epsi, g, mds);
			NEWmdVthr_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
			if (gaussian)
				NEWmdGaussian_nr(j, mu, k, xLoop, p.Epsi, g, mds);
			if (disjoint) {
				NEWmdVthrDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
				NEWmdVthrDisjointLR_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
				if (gaussian) {
					NEWmdGaussianThermal_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
					NEWmdGaussianThermal2_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
					NEWmdGaussianDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
					NEWmdGaussianLRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
					NEWmdGaussianThermalDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
					NEWmdGaussianThermalLRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
					NEWmdGaussianThermal2Disjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
					NEWmdGaussianThermal2LRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, g, mds);
				}
			}
			for (uint nu=0; nu<dim; nu++) {
				NEWddVor_nr(j, mu, k, nu, xLoop, p.Epsi, g, dds);
				NEWddVthr_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
				if (gaussian)
					NEWddGaussian_nr(j, mu, k, nu, xLoop, p.Epsi, g, dds);
				if (disjoint) {
					NEWddVthrDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
					NEWddVthrDisjointLR_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
					if (gaussian) {
						NEWddGaussianThermal_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
						NEWddGaussianThermal2_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
						NEWddGaussianDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
						NEWddGaussianLRDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
						NEWddGaussianThermalDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
						NEWddGaussianThermalLRDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
						NEWddGaussianThermal2Disjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
						NEWddGaussianThermal2LRDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, g, dds);
					}
				}
			}
		}
	}
}

cout << "-------------------------------------------------------------------------------------------------------------------------" << endl;
cout << "v middle = " << v << endl;
cout << "repulsion middle = " << repulsion << endl;
cout << "norm(mds) middle: " << mds.norm() << endl;
cout << "norm(dds) middle: " << dds.norm() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;
time = clock();

// old ones
for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dim; mu++) {
		for (uint k=0; k<N; k++) {
			if (mu==0) {
				Vor(j,k,xLoop,p.Epsi,-g,v);
				//Ver(j,k,xLoop,p.Epsi,-g,v);
				//Vdr(j,k,xLoop,p.Epsi,-g,v);
				Vthr(j,k,xLoop,beta,p.Epsi,-g,v);
				Gaussian(j, k, xLoop, p.Epsi, -g, repulsion);
				if (disjoint) {
					VthrDisjoint(j,k,xLoop,beta,p.Epsi,-g,v);
					GaussianDisjoint(j, k, xLoop, beta, p.Epsi, -g, repulsion);
					GaussianLRDisjoint(j, k, xLoop, beta, p.Epsi, -g, repulsion);
					GaussianThermal(j, k, xLoop, beta, p.Epsi, -g, repulsion);
					GaussianThermal2(j, k, xLoop, beta, p.Epsi, -g, repulsion);
					GaussianThermalDisjoint(j, k, xLoop, beta, p.Epsi, -g, repulsion);
					GaussianThermal2Disjoint(j, k, xLoop, beta, p.Epsi, -g, repulsion);
					GaussianThermalLRDisjoint(j, k, xLoop, beta, p.Epsi, -g, repulsion);
					GaussianThermal2LRDisjoint(j, k, xLoop, beta, p.Epsi, -g, repulsion);
				}
			}
			mdVor_nr(j, mu, k, xLoop, p.Epsi, -g, mds);
			mdVthr_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
			if (gaussian)
				mdGaussian_nr(j, mu, k, xLoop, p.Epsi, -g, mds);
			if (disjoint) {
				mdVthrDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
				mdVthrDisjointLR_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
				if (gaussian) {
					mdGaussianThermal_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
					mdGaussianThermal2_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
					mdGaussianDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
					mdGaussianLRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
					mdGaussianThermalDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
					mdGaussianThermal2Disjoint_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
					mdGaussianThermalLRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
					mdGaussianThermal2LRDisjoint_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
				}
			}
			for (uint nu=0; nu<dim; nu++) {
				ddVor_nr(j, mu, k, nu, xLoop, p.Epsi, -g, dds);
				ddVthr_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
				if (gaussian)
					ddGaussian_nr(j, mu, k, nu, xLoop, p.Epsi, -g, dds);
				if (disjoint) {
					ddVthrDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
					ddVthrDisjointLR_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
					if (gaussian) {
						ddGaussianThermal_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
						ddGaussianThermal2_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
						ddGaussianDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
						ddGaussianLRDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
						ddGaussianThermalDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
						ddGaussianThermalLRDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
						ddGaussianThermal2Disjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
						ddGaussianThermal2LRDisjoint_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
					}
				}
			}
		}
	}
}

cout << "-------------------------------------------------------------------------------------------------------------------------" << endl;
cout << "v after = " << v << endl;
cout << "repulsion after = " << repulsion << endl;
cout << "norm(mds) after: " << mds.norm() << endl;
cout << "norm(dds) after: " << dds.norm() << endl;
cout << "mds.size() = " << mds.size() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;

cout << "before and after should be zero, middle nonzero." << endl << endl;
cout << "==========================================================================================================" << endl;

return 0;
}
