/*
	test of nrloop stuff
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include "genloop.h"
#include "simple.h"
#include "nrloop.h"
#include "print.h"

using namespace std;

int main() {

cout << "testnrloop:" << endl;

uint dim = 4, K=4;

vec v(1);
v[0] = 2.0;

saveVectorAscii("temp/v.dat",v);

vec l;
loadVectorBinary("data/s0/loops/dim_4/K_4/loop_run_0.dat",l);
cout << "l.size() = " << l.size() << endl;

Loop<4> loop(K,0);
loop.load("data/s0/loops/dim_4/K_4/loop_run_0.dat");
cout << "loop.size() = " << loop.size() << endl;

cout << "(l-loop).(l-loop) = ";
number result = 0.0;
for (uint j=0; j<pow(2,K); j++) {
	for (uint mu = 0; mu<dim; mu++) {
		result += pow((l[dim*j+mu]-(loop[j])[mu]),2);
	}
}
cout << result << endl;

return 0;
}
