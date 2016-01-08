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

cout << "testV2" << endl;

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

number a = 0.1;
number v1r = V1r(loop,a);
number v2r = V2r(loop,a);
cout << "V1r = " << v1r << endl;
cout << "V2r = " << v2r << endl;
cout << "V2r-V1r = " << v2r-v1r << endl << endl;

number c = V1r(loop_new,a)-v1r;
number d = DV1r(loop,p,loc,a);

cout << "V1r'-V1r = " << c << endl;
cout << "DV1r = " << d << endl;
cout << "difference = " << c-d << endl << endl;

c = V2r(loop_new,a)-v2r;
d = DV2r(loop,p,loc,a);

cout << "V2r'-V2r = " << c << endl;
cout << "DV2r = " << d << endl;
cout << "difference = " << c-d << endl << endl;


return 0;
}
