/*
	testAnalysis2
		quick test of analysis.cc, the newton raphson part
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
#include "parameters.h"
#include "simple.h"
#include "print.h"

using namespace std;

int main() {

cout << "test of newton raphson part of analysis:" << endl;

cout << "-------------------------------------------------------------------------------------------" << endl;
cout << "testing NewtonRaphsonDatum:" << endl;

vector<string> ida(3);
ida[0] = "012345678912";
ida[1] = "1";
ida[2] = "10";

Parameters pa;
pa.load("inputs4");

vector<number> dataa(2);
dataa[0] = 0.123;
dataa[1] = 0.234;

cout << "asigning datum" << endl;
NewtonRaphsonDatum nrda(ida,pa,dataa), nrdb(3,2), nrdc(3,2);
string file = "data/temp/nrda.txt";
cout << "saving datum to " << file << endl;
nrda.save(file);
cout << "nrda = " << nrda << endl;

cout << "loading datum" << endl;
nrdb.load(file);

cout << "nrdb = " << nrdb << endl;

cout << "checkID = " << (int)nrdb.checkID(ida) << endl;
cout << "checkParameters = " << (int)nrdb.checkParameters(pa) << endl;
cout << "(nrda==nrdb) = " << (int)(nrda==nrdb) << endl;

vector<string> idb(3);
idb[0] = "112345678912";
idb[1] = "1";
idb[2] = "10";

nrdb.load(idb,pa,dataa);

cout << "checkID = " << (int)nrdb.checkID(ida) << endl;
cout << "checkParameters = " << (int)nrdb.checkParameters(pa) << endl;
cout << "(nrda==nrdb) = " << (int)(nrda==nrdb) << endl;

ida[0] = "";
cout << "with ida[0] empty" << endl;
cout << "checkID = " << (int)nrdb.checkID(ida) << endl;

ifstream is;
is.open(file.c_str());
is >> nrdc;
is.close();
cout << "nrdc = " << nrdc << endl;

cout << "-------------------------------------------------------------------------------------------" << endl;
cout << "testing NewtonRaphsonData:" << endl;

cout << "loading data from " << file << endl;
NewtonRaphsonData Nrd(file,3,2);

cout << "size = " << Nrd.size() << endl;

string file2 = "data/temp/nrd.txt";
cout << "saving data to " << file2 << endl;
Nrd.save(file2);

cout << "loading data from " << file2 << endl;
NewtonRaphsonData Nrd2(file2,3,2);
cout << "size = " << Nrd2.size() << endl;

cout << "find(ida) = " << (int)Nrd.find(ida) << endl;
cout << "find(pa) = " << (int)Nrd.find(pa) << endl;
cout << "find(nrdb) = " << (int)Nrd.find(nrdb) << endl;
cout << "find(ida,pa) = " << (int)Nrd.find(ida,pa) << endl;

vector<string> idc(3);
idc[0] = "";
idc[1] = "1";
idc[2] = "";

cout << "find(idc,pa) = " << (int)Nrd.find(idc,pa) << endl;

return 0;
}
