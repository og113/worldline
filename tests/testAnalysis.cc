/*
	testAnalysis
		quick test of analysis.cc
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <gsl/gsl_sf_exp.h>
#include "analysis.h"
#include "simple.h"
#include "print.h"

using namespace std;

int main() {

vector<number> v, w;

uint Size = 1000;

v.resize(Size);
number mean, mean_mcd;
number meanSqrd, meanSqrd_mcd;
number intCorrTime, expCorrTime, corrErrorSqrd, bootstrap1, bootstrap2, bootstrap3;
uint Bootstraps = 1;

srand (time(NULL));

v[0] = 5.0 + randDouble(-1.0,1.0);
mean = v[0];
meanSqrd = v[0]*v[0];
for (uint k=1; k<Size; k++) {
	v[k] = v[k-1] + ((v[k-1]-5.0)*randDouble(-0.5,0) + randDouble(-0.2,0.2));
	mean += v[k];
	meanSqrd += v[k]*v[k];
}
mean /= (number)v.size();
meanSqrd /= (number)v.size();

cout << "size = " << v.size() << endl;
cout << "mean = " << mean << endl;
cout << "meanSqrd = " << meanSqrd << endl;

saveVectorBinary<number>("data/temp/temp1.dat",v);

MonteCarloData mcd("data/temp/temp1.dat");

mcd.calcMeans(mean_mcd,meanSqrd_mcd);
mcd.calcCorrs(w,intCorrTime, expCorrTime, corrErrorSqrd);
mcd.saveResults("data/temp/temp2.dat");

ofstream ov("data/temp/v1.dat");
for (uint k=0; k<v.size(); k++)
	ov << v[k] << endl;
ov.close();

ofstream osv("data/temp/var1.dat");
for (uint k=0; k<v.size(); k++)
	osv << v[k]*v[k]-mean*mean << endl;
osv.close();

ofstream os("data/temp/corr1.dat");
for (uint k=0; k<w.size(); k++)
	os << w[k] << endl;
os.close();

uint Seed = time(NULL)+2;
bootstrap1 = mcd.calcBootstrap(Bootstraps,Seed);
bootstrap2 = mcd.calcBootstrap(2*Bootstraps,Seed);
bootstrap3 = mcd.calcBootstrap(1000*Bootstraps,Seed);

cout << "(v[v.size()-1]-5.0)/5.0 = " << (v[v.size()-1]-5.0)/5.0 << endl;
cout << "mean = " << mean << ", " << mean_mcd << endl;
cout << "meanSqrd = " << meanSqrd << ", " << meanSqrd_mcd << endl;
cout << "intCorrTime = " << intCorrTime << endl;
cout << "expCorrTime = " << expCorrTime << endl;
cout << "corrErrorSqrd = " << corrErrorSqrd << endl;
cout << "boostrapError(" << Bootstraps << ") = " << sqrt(bootstrap1) << endl;
cout << "boostrapError(" << 2*Bootstraps << ") = " << sqrt(bootstrap2) << endl;
cout << "boostrapError(" << 1000*Bootstraps << ") = " << sqrt(bootstrap3) << endl;

return 0;
}
