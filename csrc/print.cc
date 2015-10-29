/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
 
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "folder.h"
#include "simple.h" // for countLines
#include "print.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. save
	2. load
	3. explicit instantiation
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. save
-------------------------------------------------------------------------------------------------------------------------*/

// save - saveVectorAscii
template <class T>
void saveVectorAscii(const string& f,  const vector<T>& v) {
	ofstream os;
	os.open(f.c_str());
	if (os.good()) {
		os << setprecision(16);
		for (uint j=0; j<v.size(); j++) {
			os << v[j] << endl;
		}
		os.close();
	}
	else {
		cerr << "saveVectorAscii error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// save - saveVectorAsciiAppend
template <class T>
void saveVectorAsciiAppend(const string& f,  const vector<T>& v) {
	ofstream os;
	os.open(f.c_str(), ios::app);
	if (os.good()) {
		os << setprecision(16);
		for (uint j=0; j<v.size(); j++) {
			os << v[j] << endl;
		}
		os.close();
	}
	else {
		cerr << "saveVectorAsciiAppend error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// save - saveVectorBinary
template <class T>
void saveVectorBinary(const string& f,  const vector<T>& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	const T* r;
	if (os.good()) {
		for (uint j=0; j<v.size(); j++) {
			r = &v[j];
			os.write(reinterpret_cast<const char*>(r),sizeof(T));
		}
		os.close();
	}
	else {
		cerr << "saveVectorBinary error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// save - saveVectorBinaryAppend
template <class T>
void saveVectorBinaryAppend(const string& f,  const vector<T>& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary | ios::app);
	const T* r;
	if (os.good()) {
		for (uint j=0; j<v.size(); j++) {
			r = &v[j];
			os.write(reinterpret_cast<const char*>(r),sizeof(T));
		}
		os.close();
	}
	else {
		cerr << "saveVectorBinaryAppend error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. load
-------------------------------------------------------------------------------------------------------------------------*/

// loadVectorAscii
template <class T>
void loadVectorAscii(const string& f, vector<T>& v) {
	uint lines = countLines(f);
	v.resize(lines);
	ifstream is;
	is.open(f.c_str());
	if (is.good()) {
		for (uint j=0; j<lines; j++) {
			is >> v[j];
		}
		is.close();
	}
	else {
		cerr << "loadVectorAscii error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

// loadVectorBinary
template <class T>
void loadVectorBinary(const string& f, vector<T>& v) {
	T t;
	uint lines = countType<T>(f, t);
	v.resize(lines);
	ifstream is;
	is.open(f.c_str(),ios::binary);
	if (is.good()) {
		for (uint j=0; j<lines; j++) {
			is.read(reinterpret_cast<char*>(&t),sizeof(t));
			v[j] = t;
		}
		is.close();
	}
	else {
		cerr << "loadVectorBinary error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

// loadVectorAsciiColumn
template <class T>
void loadVectorAsciiColumn(const string& f, vector<T>& v, const uint& col) {
	uint cols = countColumns(f);
	if (col>cols) {
		cerr << "loadVectorAsciiColumn error: col(" << col << ")>cols(" << cols << ")" << endl;
		return;
	}
	uint lines = countLines(f);
	v.resize(lines);
	ifstream is;
	is.open(f.c_str());
	if (is.good()) {
		T dross;
		string line;
		for (uint j=0; j<lines; j++) {
			getline(is,line);
			stringstream ss(line);
			for (uint k=0; k<(col-1); k++)
				ss >> dross;
			ss >> v[j];
		}
		is.close();
	}
	else {
		cerr << "loadVectorAsciiColumn error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. explicit instantiation
-------------------------------------------------------------------------------------------------------------------------*/

// save
template void saveVectorBinary<number>(const string& f, const vector<number>& v);
template void saveVectorBinaryAppend<number>(const string& f, const vector<number>& v);
template void saveVectorAscii<number>(const string& f, const vector<number>& v);
template void saveVectorAsciiAppend<number>(const string& f, const vector<number>& v);

// load
template void loadVectorBinary<number>(const string& f, vector<number>& v);
template void loadVectorAscii<number>(const string& f, vector<number>& v);
template void loadVectorAsciiColumn<number>(const string& f, vector<number>& v, const uint& col);

