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

cout << endl << "v = ";
for (uint k=0; k<v.size(); k++)
	cout << v[k] << " ";
cout << endl;

cout << endl << "w = ";
for (uint k=0; k<w.size(); k++)
	cout << w[k] << " ";
cout << endl;

vector<string> s(3), t;
s[0] = "a    ";
s[1] = "  silly";
s[2] = "billy";

string stringFile = "data/temp/stringFile.txt";

saveVectorCsvAppend(stringFile,s);

loadVectorCsvAppend(stringFile,t);

cout << "test save and load vector csv append:" << endl;
cout << "s = ";
for (uint j=0; j<s.size(); j++)
	cout << s[j] << " ";
cout << endl << "t = ";
for (uint j=0; j<t.size(); j++)
	cout << t[j] << " ";
cout << endl;

return 0;
}
