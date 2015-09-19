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
	2 - ParametersRange
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. Parameters member functions
		- operator<<
		- save
		- load
		- empty
		- operator==
-------------------------------------------------------------------------------------------------------------------------*/

// operator<<
ostream& operator<<(ostream& os, const Parameters& p) {
	os << left;
	os << setw(20) << "LoopNum" << setw(20) << p.LoopNum << endl;
	os << setw(20) << "Ng" << setw(20) << p.Ng << endl;
	os << setw(20) << "Nms" << setw(20) << p.Nms << endl;
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
		cerr << "Parameters::load error:" << endl;
		cerr << "stream for " << filename << " not good" << endl;
		return;
	}
	string dross;
	is >> dross >> LoopNum;
	is >> dross >> Ng;
	is >> dross >> Nms;
	is >> dross >> K;
	is >> dross >> g;
	is.close();
}

// empty
bool Parameters::empty() const {
	return (LoopNum==0 && Ng==0 && Nms==0 && K==0 && abs(g)<MIN_NUMBER);
}

// operator==
bool operator==(const Parameters& l, const Parameters& r){
	return (l.LoopNum==r.LoopNum && l.Ng==r.Ng && l.Nms==r.Nms && l.K==r.K && abs(l.g-r.g)<MIN_NUMBER);
}

// writeBinary
ostream& Parameters::writeBinary(ostream& os) const {
	os.write(reinterpret_cast<const char*>(&LoopNum),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Ng),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nms),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&K),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&g),sizeof(number));
	return os;
}

// readBinary
istream& Parameters::readBinary(istream& is) {
	is.read(reinterpret_cast<char*>(&LoopNum),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Ng),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nms),sizeof(uint));
	is.read(reinterpret_cast<char*>(&K),sizeof(uint));
	is.read(reinterpret_cast<char*>(&g),sizeof(number));
	return is;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. ParametersRange
		- ParametersRange
		- operator<<
		- operator==
	
	N.B. for structs the compiler automatically writes the following member functions, unless user provided:
	 empty constructor; copy constructor; assignment operator=; destructor
-------------------------------------------------------------------------------------------------------------------------*/

// save
void ParametersRange::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		cerr << "ParametersRange::save error:" << endl;
		cerr << "stream for " << filename << " not good." << endl;
		return;
	}
	os << *this;
	os << endl;
	os.close();
}

// load
void ParametersRange::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		cerr << "ParametersRange::load error:" << endl;
		cerr << "stream for " << filename << " not good" << endl;
		return;
	}
	string dross;
	// THIS NEEDS REDOING SO THAT THERE AREN'T MISTAKES WHEN SPACES DISAPPEAR
	is >> dross >> dross >> Min.LoopNum >> dross >> Step.LoopNum >> dross >> Max.LoopNum >> dross;
	is >> dross >> dross >> Min.Ng >> dross >> Step.Ng >> dross >> Max.Ng >> dross;
	is >> dross >> dross >> Min.Nms >> dross >> Step.Nms >> dross >> Max.Nms >> dross;
	is >> dross >> dross >> Min.K >> dross >> Step.K >> dross >> Max.K >> dross;
	is >> dross >> dross >> Min.g >> dross >> Step.g >> dross >> Max.g >> dross;
	is.close();
}

// empty
bool ParametersRange::empty() const {
	return (Min.empty() && Step.empty() && Max.empty());
}

// writeBinary
ostream& ParametersRange::writeBinary(ostream& os) const {
	Min.writeBinary(os);
	Step.writeBinary(os);
	Max.writeBinary(os);
	return os;
}

// readBinary
istream& ParametersRange::readBinary(istream& is) {
	Min.readBinary(is);
	Step.readBinary(is);
	Max.readBinary(is);
	return is;
}

// operator<<
ostream& operator<<(ostream& os, const ParametersRange& pr) {
	os << left;
	os << setw(20) << "LoopNum" << setw(20) << "[ " << (pr.Min).LoopNum << " , " \
						<< (pr.Step).LoopNum << " , " << (pr.Max).LoopNum << " ]" << endl;
	os << setw(20) << "Ng" << setw(20) << "[ " << (pr.Min).Ng << " , " \
						<< (pr.Step).Ng << " , " << (pr.Max).Ng << " ]" << endl;
	os << setw(20) << "Nms" << setw(20) << "[ " << (pr.Min).Nms << " , " \
						<< (pr.Step).Nms << " , " << (pr.Max).Nms << " ]" << endl;
	os << setw(20) << "K" << setw(20) << "[ " << (pr.Min).K << " , " \
						<< (pr.Step).K << " , " << (pr.Max).K << " ]" << endl;
	os << setw(20) << "g" << setw(20) << "[ " << (pr.Min).g << " , " \
						<< (pr.Step).g << " , " << (pr.Max).g << " ]" << endl;
	return os;
}

// operator==
bool operator==(const ParametersRange& l, const ParametersRange& r) {
	return (l.Min==r.Min && l.Step==r.Step && l.Max==r.Step);
}

