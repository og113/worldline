/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "simple.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - Parameters
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	2. Parameters member functions
		- operator<<
		- save
		- load
		- empty
		- operator==
-------------------------------------------------------------------------------------------------------------------------*/

// operator<<
ostream& operator<<(ostream& os, const Parameters& p) {
	os << left;
	os << setw(20) << "Loops" << setw(20) << p.Loops << endl;
	os << setw(20) << "K" << setw(20) << p.K << endl;
	os << setw(20) << "g" << setw(20) << p.g << endl;
	return os;
}

//save
void Parameters::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		cerr << "Parameters::save error:" << endl;
		cerr << "stream for " << filename << " not good." << endl;
		return;
	}
	os << *this;
	os << endl;
	os.close();
}

//load
void Parameters::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		cerr << "PrimaryParameters::load error:" << endl;
		cerr << "stream for " << filename << " not good" << endl;
		return;
	}
	string dross;
	is >> dross >> Loops;
	is >> dross >> K;
	is >> dross >> g;
	is.close();
}

// empty
bool Parameters::empty() const {
	return (Loops==0 && K==0 && abs(g)<MIN_NUMBER);
}

// operator==
bool operator==(const Parameters& l, const Parameters& r){
	return (l.Loops==r.Loops && l.K==r.K && abs(l.g-r.g)<MIN_NUMBER);
}

// writeBinary
ostream& Parameters::writeBinary(ostream& os) const {
	os.write(reinterpret_cast<const char*>(&Loops),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&K),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&g),sizeof(number));
	return os;
}

// readBinary
istream& Parameters::readBinary(istream& is) {
	is.read(reinterpret_cast<char*>(&Loops),sizeof(uint));
	is.read(reinterpret_cast<char*>(&K),sizeof(uint));
	is.read(reinterpret_cast<char*>(&g),sizeof(number));
	return is;
}

