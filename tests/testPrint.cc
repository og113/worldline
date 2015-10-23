/*
	testPrint
		quick test of print.cc
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "simple.h"
#include "print.h"

using namespace std;

int main() {

vector<number> v, w;

v.resize(8,6.6);

saveVectorBinary<number>("data/temp/temp1.dat",v);

loadVectorBinary<number>("data/temp/temp1.dat",w);

for (uint k=0; k<w.size(); k++)
	cout << w[k] << endl;

return 0;
}
