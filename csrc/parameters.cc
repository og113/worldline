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

// Size
const uint Parameters::Size = 5;
const uint ParametersRange::Size = 5;

// step
void Parameters::step(const ParametersRange& pr) {
	Label label = static_cast<Parameters::Label>(0);
	if (pr.toStep(label)) {
		number stepSize;
		switch (label){
			case loops:
				stepSize = ((pr.Max).Loops-(pr.Min).Loops)/((pr.Steps)[label]-1.0);
				Loops += (uint)stepSize;
				break;
			case ng:
				stepSize = ((pr.Max).Ng-(pr.Min).Ng)/((pr.Steps)[label]-1.0);
				Ng += (uint)stepSize;
				break;
			case nms:
				stepSize = ((pr.Max).Nms-(pr.Min).Nms)/((pr.Steps)[label]-1.0);
				Nms += (uint)stepSize;
				break;
			case k:
				stepSize = ((pr.Max).K-(pr.Min).K)/((pr.Steps)[label]-1.0);
				K += (uint)stepSize;
				break;
			case g:
				stepSize = ((pr.Max).G-(pr.Min).G)/((pr.Steps)[label]-1.0);
				G += stepSize;
				break;
			default:
				cerr << "Parameters error: Label unknown" << endl;
				break;
		}
	}
	else {
		cerr << "Parameters error: cannot step as Steps vector does not contain one nonZero element" << endl;
	}
}

// operator<<
ostream& operator<<(ostream& os, const Parameters& p) {
	os << left;
	os << setw(20) << "Loops" << setw(20) << p.Loops << endl;
	os << setw(20) << "Ng" << setw(20) << p.Ng << endl;
	os << setw(20) << "Nms" << setw(20) << p.Nms << endl;
	os << setw(20) << "K" << setw(20) << p.K << endl;
	os << setw(20) << "G" << setw(20) << p.G << endl;
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
	is >> dross >> Loops;
	is >> dross >> Ng;
	is >> dross >> Nms;
	is >> dross >> K;
	is >> dross >> G;
	is.close();
}

// empty
bool Parameters::empty() const {
	return (Loops==0 && Ng==0 && Nms==0 && K==0 && abs(g)<MIN_NUMBER);
}

// operator==
bool operator==(const Parameters& l, const Parameters& r){
	return (l.Loops==r.Loops && l.Ng==r.Ng && l.Nms==r.Nms && l.K==r.K && abs(l.G-r.G)<MIN_NUMBER);
}

// writeBinary
ostream& Parameters::writeBinary(ostream& os) const {
	os.write(reinterpret_cast<const char*>(&Loops),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Ng),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nms),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&K),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&G),sizeof(number));
	return os;
}

// readBinary
istream& Parameters::readBinary(istream& is) {
	is.read(reinterpret_cast<char*>(&Loops),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Ng),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nms),sizeof(uint));
	is.read(reinterpret_cast<char*>(&K),sizeof(uint));
	is.read(reinterpret_cast<char*>(&G),sizeof(number));
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

// ParametersRange empty constructor
ParametersRange::ParametersRange(): Min(), Max() {
	Steps = vector<uint>(Size,0);
}

// ParametersRange constructor
ParametersRange::ParametersRange(const Parameters& min, const Parameters& max, const vector<uint>& steps): \
									Min(min), Max(max) {
	if (steps.size()==Size)
		Steps = steps;
	else {
		Steps = vector<uint>(Size,0);
		cerr << "ParametersRange Error: Initialization not possible as input steps.size() = " << steps.size() << " not " << Size << endl;
	} 							
}

// toStep
bool ParametersRange::toStep(Parameters::Label& stepNum) const {
	uint nonZeros = 0;
	for (uint j=0; j<Size; j++) {
		if (Steps[j]>0) {
			stepNum = static_cast<Parameters::Label>(j);
			nonZeros++;
		}
	}
	return (nonZeros==1);
}

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
	is >> dross >> dross >> Min.Loops >> dross >> Max.Loops >> dross >> Steps[0] >> dross;
	is >> dross >> dross >> Min.Ng >> dross >> Max.Ng >> dross >> Steps[1] >> dross;
	is >> dross >> dross >> Min.Nms >> dross >> Max.Nms >> dross >> Steps[2] >> dross;
	is >> dross >> dross >> Min.K >> dross >> Max.K >> dross >> Steps[3] >> dross;
	is >> dross >> dross >> Min.G >> dross >> Max.G >> dross >> Steps[4] >> dross;
	is.close();
}

// empty
bool ParametersRange::empty() const {
	return (Min.empty() && Max.empty() && Steps.empty());
}

// writeBinary
ostream& ParametersRange::writeBinary(ostream& os) const {
	Min.writeBinary(os);
	Max.writeBinary(os);
	os.write(reinterpret_cast<const char*>(&Steps[0]),Size*sizeof(uint));
	return os;
}

// readBinary
istream& ParametersRange::readBinary(istream& is) {
	Min.readBinary(is);
	Max.readBinary(is);
	is.read(reinterpret_cast<char*>(&Steps[0]),Size*sizeof(uint));
	return is;
}

// operator<<
ostream& operator<<(ostream& os, const ParametersRange& pr) {
	os << left;
	os << setw(20) << "Loops" << "[ " << setw(12) << (pr.Min).Loops << " , " \
						<< setw(12) << (pr.Max).Loops << " , " << setw(12) << (pr.Steps)[0] << " ]" << endl;
	os << setw(20) << "Ng" << "[ " << setw(12) <<  (pr.Min).Ng << " , " \
						<< setw(12) << (pr.Max).Ng << " , " << setw(12) << (pr.Steps)[1] << " ]" << endl;
	os << setw(20) << "Nms" << "[ " << setw(12) << (pr.Min).Nms << " , " \
						<< setw(12) << (pr.Max).Nms << " , " << setw(12) << (pr.Steps)[2] << " ]" << endl;
	os << setw(20) << "K" << "[ " << setw(12) << setw(12) << (pr.Min).K << " , " \
						<< setw(12) << (pr.Max).K << " , " << setw(12) << (pr.Steps)[3] << " ]" << endl;
	os << setw(20) << "g" << "[ " << setw(12) << (pr.Min).G << " , " \
						<< setw(12) << (pr.Max).G << " , " << setw(12) << (pr.Steps)[4] << " ]" << endl;
	return os;
}

// operator==
bool operator==(const ParametersRange& l, const ParametersRange& r) {
	return (l.Min==r.Min && l.Max==r.Max && l.Steps==r.Steps);
}

