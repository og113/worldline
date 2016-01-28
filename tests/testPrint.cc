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

number Length(const vector<number>& v) {
	number result = 0.0;
	for (uint k=0; k<v.size(); k++) {
		result += v[k]*v[k];
	}
	return sqrt(result);
}

int main() {

vector<number> v, w;

v.resize(8,6.6);

cout << "Test Print passes if following numbers agree: " << endl;

saveVectorBinary("data/temp/temp1.dat",v);

loadVectorBinary("data/temp/temp1.dat",w);
cout << v.size() << " " << w.size() << endl;
cout << Length(v) << " " << Length(w) << endl;

cout << endl << "v = " << endl;
for (uint k=0; k<v.size(); k++)
	cout << v[k] << endl;

cout << endl << "w = " << endl;
for (uint k=0; k<w.size(); k++)
	cout << w[k] << endl;

return 0;
}
