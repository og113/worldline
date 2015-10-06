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
const uint Parameters::Size = 8;
const uint ParametersRange::Size = 8;

// step
void Parameters::step(const ParametersRange& pr) {
	Label label = static_cast<Parameters::Label>(0);
	if (pr.toStep(label)) {
		number stepSize;
		switch (label){
			case nl:
				stepSize = ((pr.Max).Nl-(pr.Min).Nl)/((pr.Steps)[label-1]-1.0);
				Nl += (uint)stepSize;
				break;
			case ng:
				stepSize = ((pr.Max).Ng-(pr.Min).Ng)/((pr.Steps)[label-1]-1.0);
				Ng += (uint)stepSize;
				break;
				break;
			case nig:
				stepSize = ((pr.Max).Nig-(pr.Min).Nig)/((pr.Steps)[label-1]-1.0);
				Nig += (uint)stepSize;
				break;
			case nsw:
				stepSize = ((pr.Max).Nsw-(pr.Min).Nsw)/((pr.Steps)[label-1]-1.0);
				Nsw += (uint)stepSize;
			case npsw:
				stepSize = ((pr.Max).Npsw-(pr.Min).Npsw)/((pr.Steps)[label-1]-1.0);
				Npsw += (uint)stepSize;
				break;
			case k:
				stepSize = ((pr.Max).K-(pr.Min).K)/((pr.Steps)[label-1]-1.0);
				K += (uint)stepSize;
				break;
			case g:
				stepSize = ((pr.Max).G-(pr.Min).G)/((pr.Steps)[label-1]-1.0);
				G += stepSize;
				break;
			case epsi:
				stepSize = ((pr.Max).Epsi-(pr.Min).Epsi)/((pr.Steps)[label-1]-1.0);
				Epsi += stepSize;
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
	os << setw(20) << "Nl" << setw(20) << p.Nl << endl;
	os << setw(20) << "Ng" << setw(20) << p.Ng << endl;
	os << setw(20) << "Nig" << setw(20) << p.Nig << endl;
	os << setw(20) << "Nsw" << setw(20) << p.Nsw << endl;
	os << setw(20) << "Npsw" << setw(20) << p.Npsw << endl;
	os << setw(20) << "K" << setw(20) << p.K << endl;
	os << setw(20) << "G" << setw(20) << p.G << endl;
	os << setw(20) << "Epsi" << setw(20) << p.Epsi << endl;
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
	is >> dross >> Nl;
	is >> dross >> Ng;
	is >> dross >> Nig;
	is >> dross >> Nsw;
	is >> dross >> Npsw;
	is >> dross >> K;
	is >> dross >> G;
	is >> dross >> Epsi;
	is.close();
}

// empty
bool Parameters::empty() const {
	return (Nl==0 && Ng==0 && Nig==0 && Nsw==0 && Npsw==0 && K==0 && abs(G)<MIN_NUMBER);
}

// operator==
bool operator==(const Parameters& l, const Parameters& r){
	return (l.Nl==r.Nl && l.Ng==r.Ng && l.Nig==r.Nig && && l.Nsw==r.Nsw && && l.Npsw==r.Npsw && l.K==r.K && abs(l.G-r.G)<MIN_NUMBER);
}

// writeBinary
ostream& Parameters::writeBinary(ostream& os) const {
	os.write(reinterpret_cast<const char*>(&Nl),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Ng),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nig),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nsw),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Npsw),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&K),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&G),sizeof(number));
	os.write(reinterpret_cast<const char*>(&Epsi),sizeof(number));
	return os;
}

// readBinary
istream& Parameters::readBinary(istream& is) {
	is.read(reinterpret_cast<char*>(&Nl),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Ng),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nig),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nsw),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Npsw),sizeof(uint));
	is.read(reinterpret_cast<char*>(&K),sizeof(uint));
	is.read(reinterpret_cast<char*>(&G),sizeof(number));
	is.read(reinterpret_cast<char*>(&Epsi),sizeof(number));
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
			stepNum = static_cast<Parameters::Label>(j+1);
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
	is >> dross >> dross >> Min.Nl >> dross >> Max.Nl >> dross >> Steps[0] >> dross;
	is >> dross >> dross >> Min.Ng >> dross >> Max.Ng >> dross >> Steps[1] >> dross;
	is >> dross >> dross >> Min.Nig >> dross >> Max.Nig >> dross >> Steps[3] >> dross;
	is >> dross >> dross >> Min.Nsw >> dross >> Max.Nsw >> dross >> Steps[2] >> dross;
	is >> dross >> dross >> Min.Npsw >> dross >> Max.Npsw >> dross >> Steps[4] >> dross;
	is >> dross >> dross >> Min.K >> dross >> Max.K >> dross >> Steps[5] >> dross;
	is >> dross >> dross >> Min.G >> dross >> Max.G >> dross >> Steps[6] >> dross;
	is >> dross >> dross >> Min.Epsi >> dross >> Max.Epsi >> dross >> Steps[7] >> dross;
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
	os << setw(20) << "Nl" << "[ " << setw(12) << (pr.Min).Nl << " , " \
						<< setw(12) << (pr.Max).Nl << " , " << setw(12) << (pr.Steps)[0] << " ]" << endl;
	os << setw(20) << "Ng" << "[ " << setw(12) <<  (pr.Min).Ng << " , " \
						<< setw(12) << (pr.Max).Ng << " , " << setw(12) << (pr.Steps)[1] << " ]" << endl;
	os << setw(20) << "Nig" << "[ " << setw(12) << (pr.Min).Nig << " , " \
						<< setw(12) << (pr.Max).Nig << " , " << setw(12) << (pr.Steps)[3] << " ]" << endl;
	os << setw(20) << "Nsw" << "[ " << setw(12) << (pr.Min).Nsw << " , " \
						<< setw(12) << (pr.Max).Nsw << " , " << setw(12) << (pr.Steps)[2] << " ]" << endl;
	os << setw(20) << "Npsw" << "[ " << setw(12) << (pr.Min).Npsw << " , " \
						<< setw(12) << (pr.Max).Npsw << " , " << setw(12) << (pr.Steps)[4] << " ]" << endl;
	os << setw(20) << "K" << "[ " << setw(12) << setw(12) << (pr.Min).K << " , " \
						<< setw(12) << (pr.Max).K << " , " << setw(12) << (pr.Steps)[5] << " ]" << endl;
	os << setw(20) << "G" << "[ " << setw(12) << (pr.Min).G << " , " \
						<< setw(12) << (pr.Max).G << " , " << setw(12) << (pr.Steps)[6] << " ]" << endl;
	os << setw(20) << "Epsi" << "[ " << setw(12) << (pr.Min).Epsi << " , " \
						<< setw(12) << (pr.Max).Epsi << " , " << setw(12) << (pr.Steps)[7] << " ]" << endl;
	return os;
}

// operator==
bool operator==(const ParametersRange& l, const ParametersRange& r) {
	return (l.Min==r.Min && l.Max==r.Max && l.Steps==r.Steps);
}

