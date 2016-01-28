/*
	testDV
*/

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "simple.h"
#include "print.h"
#include "genloop.h"

using namespace std;

int main() {

cout << "testDV" << endl;

#define dim 4

uint K=10, Length = pow(2,K);
uint Seed = 1;
//string file = "data/s0/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(K)+"/loop_run_0.dat";
string file = "data/circle/loops/dim_"+nts<uint>(dim)+"/K_"+nts<uint>(K)+"/loop_R_1_rank_0.dat";

cout << "loading from " << file << endl;

Loop<dim> loop(K,Seed), loop_new(K,Seed);
loop.load(file);
loop_new.load(file);

uint loc = Length/2;
Point<dim> p = loop[loc];
for (uint n=0; n<dim; n++)
	p[n] += 0.2*(number)n;
	
cout << "changing point: " << endl << loop[loc] << endl;
cout << "to point: "<< endl << p << endl << endl;
	
loop_new[loc] = p;

number v1 = V1(loop);
number c = V1(loop_new)-v1;
number d = DV1(loop,p,loc);
cout << "V1 = " << v1 << endl;
cout << "V1'-V1 = " << c << endl;
cout << "DV1 = " << d << endl;
cout << "difference = " << c-d << endl << endl;

number a = 0.1;
number v1r = V1r(loop,a);
c = V1r(loop_new,a)-v1r;
d = DV1r(loop,p,loc,a);
cout << "V1r = " << v1r << endl;
cout << "V1r'-V1r = " << c << endl;
cout << "DV1r = " << d << endl;
cout << "difference = " << c-d << endl << endl;

number v0 = V0(loop);
c = V0(loop_new)-v0;
d = DV0(loop,p,loc);
cout << "V0 = " << v0 << endl;
cout << "V0'-V0 = " << c << endl;
cout << "DV0 = " << d << endl;
cout << "aprxDV0 = " << aprxDV0(loop,p,loc) << endl;
cout << "difference = " << c-d << endl << endl;

c = Sm(loop_new)-Sm(loop);
d = DSm(loop,p,loc);
cout << "Sm = " << Sm(loop) << endl;
cout << "Sm'-Sm = " << c << endl;
cout << "DSm = " << d << endl;
cout << "difference = " << c-d << endl << endl;

c = I0(loop_new)-I0(loop);
d = DI0(loop,p,loc);
cout << "I0 = " << I0(loop) << endl;
cout << "I0'-I0 = " << c << endl;
cout << "DI0 = " << d << endl;
cout << "difference = " << c-d << endl << endl;

c = FGamma(loop_new)-FGamma(loop);
d = DFGamma(loop,p,loc);
number gamma = 2.0*PI/(number)loop.size();
number temp, cot_gamma;
temp = Dot(loop[2],loop[1],loop[1],loop[0]);
cot_gamma = DistanceSquared(loop[2],loop[1])*DistanceSquared(loop[1],loop[0]) - temp*temp;
cot_gamma = temp/sqrt(cot_gamma);

cout << "gamma for circle = " << gamma << ", " << atan(1.0/cot_gamma) << endl;
cout << "FGamma per site = " << (gamma/tan(gamma)-1.0) << ", " << cot_gamma*atan(1.0/cot_gamma)-1.0 << endl;
cout << "FGamma for circle = " << loop.size()*(gamma/tan(gamma)-1.0) << ", ";
cout << loop.size()*(cot_gamma*atan(1.0/cot_gamma)-1.0) << endl << endl;

cout << "FGamma = " << FGamma(loop) << endl;
cout << "FGamma'-FGamma = " << c << endl;
cout << "DFGamma = " << d << endl;
cout << "difference = " << c-d << endl << endl;

c = S0(loop_new)-S0(loop);
d = DS0(loop,p,loc);
cout << "S0 = " << S0(loop) << endl;
cout << "S0'-S0 = " << c << endl;
cout << "DS0 = " << d << endl;
cout << "difference = " << c-d << endl << endl;

cout << "DeltaGamma = " << FGamma(loop)*log(L(loop)/a) << endl;
cout << "DeltaSelfErg = " << PI*L(loop)/a << endl;
cout << "V1r - Deltas = " << v1r-PI*L(loop)/a-FGamma(loop)*log(L(loop)/a) << endl;
cout << "2*PI*PI = " << -2.0*PI*PI << endl;

return 0;
}
