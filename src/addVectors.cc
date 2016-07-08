/*
	addVectors
		quick program to load two files containing vectors of doubles, add them and print the result
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "simple.h"
#include "folder.h"
#include "genloop.h"
#include "print.h"

using namespace std;

struct FileType {
	enum Option { binary, ascii, loop};
};

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. getting argv
----------------------------------------------------------------------------------------------------------------------------*/

Filename f1, f2, fo;
string filetype = "";
uint K=0;

//dimension
#define dim 4

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("file1")==0 || id.compare("f1")==0 || id.compare("1")==0) f1 = (string)argv[2*j+2];
		else if (id.compare("file2")==0 || id.compare("f2")==0 || id.compare("2")==0) f2 = (string)argv[2*j+2];
		else if (id.compare("o")==0 || id.compare("out")==0 || id.compare("outFile")==0 || id.compare("outfile")==0)\
								 fo = (string)argv[2*j+2];
		else if (id.compare("l")==0 || id.compare("loop")==0) filetype = "loop";
		else if (id.compare("b")==0 || id.compare("binary")==0) filetype = "binary";
		else if (id.compare("a")==0 || id.compare("ascii")==0) filetype = "ascii";
		else if (id.compare("filetype")==0 || id.compare("type")==0) filetype = (string)argv[2*j+2];
		else if (id.compare("K")==0) K = stn<uint>(argv[2*j+2]);
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

FileType::Option ft = FileType::loop;
if (!filetype.empty()) {
	if (filetype.compare("binary")==0)
		ft = FileType::binary;
	else if (filetype.compare("ascii")==0)
		ft = FileType::ascii;
	else if (filetype.compare("loop")==0)
		ft = FileType::loop;
	else
		cerr << "addVectors error: file type not understood: " << filetype << endl;
}


if (f1.empty() || f2.empty()) {
	cerr << "must provide file1 and file2 as input" << endl;
	cerr << f1 << ", " << f2 << endl;
	return 1;
}

if (ft==FileType::loop && K==0) {
	vector<number> v;
	loadVectorBinary< vector<number> >(f1,v);
	uint K1=log2(v.size()/dim); // the -1 assumes a translation invariance lagrange multiplier
	loadVectorBinary< vector<number> >(f2,v);
	uint K2=log2(v.size()/dim);
	if (K1==K2)
		K = K1;
	else {
		cerr << "addVectors error: K1!=K2" << endl;
		return 1;
	}
}

if (fo.empty()) {
	fo = "data/temp/addVectors.dat";
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. loading and printing
----------------------------------------------------------------------------------------------------------------------------*/

if (ft==FileType::loop) {
	Loop<dim> l1(K,0), l2(K,0);
	l1.load(f1);
	l2.load(f2);
	for (uint j=0; j<pow(2,K); j++)
		l1[j] += l2[j];
	l1.saveAscii(fo);
	cout << "saving as loop with K=" << K << endl;
}
else if (ft==FileType::binary) {
	vector<number> v1, v2;
	loadVectorBinary< vector<number> >(f1,v1);
	loadVectorBinary< vector<number> >(f2,v2);
	if (v1.size()>pow(2,K) && K!=0)
		v1.resize(pow(2,K));
	if (v2.size()>pow(2,K) && K!=0)
		v2.resize(pow(2,K));
	for (uint j=0; j<pow(2,K); j++)
		v1[j] += v2[j];
	saveVectorAscii< vector<number> >(fo,v1);
	cout << "saving as ascii vector with K=" << K << endl;
}
else if (ft==FileType::ascii) {
	vector<number> v1, v2;
	loadVectorAscii< vector<number> >(f1,v1);
	loadVectorAscii< vector<number> >(f2,v2);
	if (v1.size()>pow(2,K) && K!=0)
		v1.resize(pow(2,K));
	if (v2.size()>pow(2,K) && K!=0)
		v2.resize(pow(2,K));
	for (uint j=0; j<pow(2,K); j++)
		v1[j] += v2[j];
	saveVectorAscii< vector<number> >(fo,v1);
	cout << "saving as ascii vector with K=" << K << endl;
}

cout << "saved sum of vectors from " << f1 << " and " << f2 << " to " << fo << " in ascii format." << endl;

return 0;
}
