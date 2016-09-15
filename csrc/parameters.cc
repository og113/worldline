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
const uint Parameters::Size = 16;
const uint ParametersRange::Size = 16;

// nameVector()
vector<string> Parameters::nameVector() const {
	vector<string> v(Size);
	v[0] = "Nl";
	v[1] =  "Ng";
	v[2] =  "Nig";
	v[3] =  "Nsw";
	v[4] =  "Npsw";
	v[5] =  "K";
	v[6] =  "G";
	v[7] =  "B";
	v[8] =  "T";
	v[9] =  "Epsi";
	v[10] =  "Mu";
	v[11] =  "P1";
	v[12] =  "P2";
	v[13] =  "P3";
	v[14] =  "P4";
	v[15] =  "Lambda";
	return v;
}

// valueVector()
vector<string> Parameters::valueVector() const {
	vector<string> v(Size);
	v[0] = nts(Nl);
	v[1] =  nts(Ng);
	v[2] =  nts(Nig);
	v[3] =  nts(Nsw);
	v[4] =  nts(Npsw);
	v[5] =  nts(K);
	v[6] =  nts(G,16);
	v[7] =  nts(B,16);
	v[8] =  nts(T,16);
	v[9] =  nts(Epsi,16);
	v[10] =  nts(Mu,16);
	v[11] =  nts(P1,16);
	v[12] =  nts(P2,16);
	v[13] =  nts(P3,16);
	v[14] =  nts(P4,16);
	v[15] =  nts(Lambda,16);
	return v;
}

// step
void Parameters::step(const ParametersRange& pr, const Parameters::Label& label, const uint& j) {
	if (pr.toStep(label)) {
		number stepSize;
		switch (label) {
			case nl:
				stepSize = ((int)((pr.Max).Nl-(pr.Min).Nl))/((pr.Steps)[label-1]-1.0);
				Nl += (int)j*stepSize;
				break;
			case ng:
				stepSize = ((int)((pr.Max).Ng-(pr.Min).Ng))/((pr.Steps)[label-1]-1.0);
				Ng += (int)j*stepSize;
				break;
			case nig:
				stepSize = ((int)((pr.Max).Nig-(pr.Min).Nig))/((pr.Steps)[label-1]-1.0);
				Nig += (int)j*stepSize;
				break;
			case nsw:
				stepSize = ((int)((pr.Max).Nsw-(pr.Min).Nsw))/((pr.Steps)[label-1]-1.0);
				Nsw += (int)j*stepSize;
			case npsw:
				stepSize = ((int)((pr.Max).Npsw-(pr.Min).Npsw))/((pr.Steps)[label-1]-1.0);
				Npsw += (int)j*stepSize;
				break;
			case k:
				stepSize = ((int)((pr.Max).K-(pr.Min).K))/((pr.Steps)[label-1]-1.0);
				K += (int)j*stepSize;
				break;
			case g:
				stepSize = ((pr.Max).G-(pr.Min).G)/((pr.Steps)[label-1]-1.0);
				G += j*stepSize;
				break;
			case b:
				stepSize = ((pr.Max).B-(pr.Min).B)/((pr.Steps)[label-1]-1.0);
				B += j*stepSize;
				break;
			case t:
				stepSize = ((pr.Max).T-(pr.Min).T)/((pr.Steps)[label-1]-1.0);
				T += j*stepSize;
				break;
			case epsi:
				stepSize = ((pr.Max).Epsi-(pr.Min).Epsi)/((pr.Steps)[label-1]-1.0);
				Epsi += j*stepSize;
				break;
			case mu:
				stepSize = ((pr.Max).Mu-(pr.Min).Mu)/((pr.Steps)[label-1]-1.0);
				Mu += j*stepSize;
				break;
			case p1:
				stepSize = ((pr.Max).P1-(pr.Min).P1)/((pr.Steps)[label-1]-1.0);
				P1 += j*stepSize;
				break;
			case p2:
				stepSize = ((pr.Max).P2-(pr.Min).P2)/((pr.Steps)[label-1]-1.0);
				P2 += j*stepSize;
				break;
			case p3:
				stepSize = ((pr.Max).P3-(pr.Min).P3)/((pr.Steps)[label-1]-1.0);
				P4 += j*stepSize;
				break;
			case p4:
				stepSize = ((pr.Max).P4-(pr.Min).P4)/((pr.Steps)[label-1]-1.0);
				P4 += j*stepSize;
				break;	
			case lambda:
				stepSize = ((pr.Max).Lambda-(pr.Min).Lambda)/((pr.Steps)[label-1]-1.0);
				Lambda += j*stepSize;
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

// step
void Parameters::step(const ParametersRange& pr, const Parameters::Label& label) {
	step(pr,label,1);
}

// operator<<
ostream& operator<<(ostream& os, const Parameters& p) {
	os << left;
	vector<string> nv = p.nameVector();
	vector<string> vv = p.valueVector();
	for (uint j=0; j<p.Size; j++)
		os << setw(20) << nv[j] << setw(20) << vv[j] << endl;
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
void Parameters::load(const vector<string>& v) {
		
	Nl = stn<uint>(v[0]);
	Ng = stn<uint>(v[1]);
	Nig = stn<uint>(v[2]);
	Nsw = stn<uint>(v[3]);
	Npsw = stn<uint>(v[4]);
	K = stn<uint>(v[5]);
	G = stn<number>(v[6]);
	B = stn<number>(v[7]);
	T = stn<number>(v[8]);
	Epsi = stn<number>(v[9]);
	Mu = stn<number>(v[10]);
	P1 = stn<number>(v[11]);
	P2 = stn<number>(v[12]);
	P3 = stn<number>(v[13]);
	P4 = stn<number>(v[14]);
	Lambda = stn<number>(v[15]);

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
	vector<string> v(Size);
	for (uint j=0; j<Size; j++)
		is >> dross >> v[j];
	is.close();
		
	load(v);

}

// empty
bool Parameters::empty() const {
	return (Nl==0 && Ng==0 && Nig==0 && Nsw==0 && Npsw==0 && K==0 && abs(G)<MIN_NUMBER && abs(B)<MIN_NUMBER\
	 && abs(T)<MIN_NUMBER && abs(Epsi)<MIN_NUMBER && abs(Mu)<MIN_NUMBER && abs(P1)<MIN_NUMBER && abs(P2)<MIN_NUMBER\
	  && abs(P3)<MIN_NUMBER && abs(P4)<MIN_NUMBER && abs(Lambda)<MIN_NUMBER);
}

// operator==
bool operator==(const Parameters& l, const Parameters& r){
	return (l.Nl==r.Nl && l.Ng==r.Ng && l.Nig==r.Nig && l.Nsw==r.Nsw && l.Npsw==r.Npsw && l.K==r.K && abs(l.G-r.G)<MIN_NUMBER\
	 && abs(l.B-r.B)<MIN_NUMBER && abs(l.T-r.T)<MIN_NUMBER && abs(l.Epsi-r.Epsi)<MIN_NUMBER && abs(l.Mu-r.Mu)<MIN_NUMBER\
	 && abs(l.P1-r.P1)<MIN_NUMBER && abs(l.P2-r.P2)<MIN_NUMBER && abs(l.P3-r.P3)<MIN_NUMBER \
	 && abs(l.Lambda-r.Lambda)<MIN_NUMBER);
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
	os.write(reinterpret_cast<const char*>(&B),sizeof(number));
	os.write(reinterpret_cast<const char*>(&T),sizeof(number));
	os.write(reinterpret_cast<const char*>(&Epsi),sizeof(number));
	os.write(reinterpret_cast<const char*>(&Mu),sizeof(number));
	os.write(reinterpret_cast<const char*>(&P1),sizeof(number));
	os.write(reinterpret_cast<const char*>(&P2),sizeof(number));
	os.write(reinterpret_cast<const char*>(&P3),sizeof(number));
	os.write(reinterpret_cast<const char*>(&P4),sizeof(number));
	os.write(reinterpret_cast<const char*>(&Lambda),sizeof(number));
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
	is.read(reinterpret_cast<char*>(&B),sizeof(number));
	is.read(reinterpret_cast<char*>(&T),sizeof(number));
	is.read(reinterpret_cast<char*>(&Epsi),sizeof(number));
	is.read(reinterpret_cast<char*>(&Mu),sizeof(number));
	is.read(reinterpret_cast<char*>(&P1),sizeof(number));
	is.read(reinterpret_cast<char*>(&P2),sizeof(number));
	is.read(reinterpret_cast<char*>(&P3),sizeof(number));
	is.read(reinterpret_cast<char*>(&P4),sizeof(number));
	is.read(reinterpret_cast<char*>(&Lambda),sizeof(number));
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
bool ParametersRange::toStep(const Parameters::Label& stepNum) const {
	return (Steps[stepNum-1]>0);
}

// totalSteps
uint ParametersRange::totalSteps() const {
	uint N = 1;
	for (uint j=0; j<ParametersRange::Size; j++)
		N *= (Steps[j]>0? Steps[j]: 1);
	return N;
}

// position
Parameters ParametersRange::position(const uint& pos) const {
	Parameters p =  Min;
	uint n, Nratio, local;
	Nratio = totalSteps();
	local = pos;
	for (uint j=0; j<ParametersRange::Size; j++) {
		if (Steps[ParametersRange::Size-j-1]>0) {
			Nratio /= Steps[ParametersRange::Size-j-1];
			n = local/Nratio;
			local -= n*Nratio;
			p.step(*this,static_cast<Parameters::Label>(ParametersRange::Size-j),n);
		}
	}
	return p;
}

// neigh
Parameters ParametersRange::neigh(const uint& pos) const {
	if (pos==0) {
		cerr << "ParametersRange::neigh error: position 0 has no neighbour" << endl;
		return Min;
	}
	uint N = 1, k=0, npos = pos;
	while (k<ParametersRange::Size && npos==pos) {
		if (Steps[k]>0) {
			if (pos%(N*Steps[k])!=0)
				npos -= N;
			N *= Steps[k];
		}
		k++;
	}
	return position(npos);
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
	is >> dross >> dross >> Min.Nig >> dross >> Max.Nig >> dross >> Steps[2] >> dross;
	is >> dross >> dross >> Min.Nsw >> dross >> Max.Nsw >> dross >> Steps[3] >> dross;
	is >> dross >> dross >> Min.Npsw >> dross >> Max.Npsw >> dross >> Steps[4] >> dross;
	is >> dross >> dross >> Min.K >> dross >> Max.K >> dross >> Steps[5] >> dross;
	is >> dross >> dross >> Min.G >> dross >> Max.G >> dross >> Steps[6] >> dross;
	is >> dross >> dross >> Min.B >> dross >> Max.B >> dross >> Steps[7] >> dross;
	is >> dross >> dross >> Min.T >> dross >> Max.T >> dross >> Steps[8] >> dross;
	is >> dross >> dross >> Min.Epsi >> dross >> Max.Epsi >> dross >> Steps[9] >> dross;
	is >> dross >> dross >> Min.Mu >> dross >> Max.Mu >> dross >> Steps[10] >> dross;
	is >> dross >> dross >> Min.P1 >> dross >> Max.P1 >> dross >> Steps[11] >> dross;
	is >> dross >> dross >> Min.P2 >> dross >> Max.P2 >> dross >> Steps[12] >> dross;
	is >> dross >> dross >> Min.P3 >> dross >> Max.P3 >> dross >> Steps[13] >> dross;
	is >> dross >> dross >> Min.P4 >> dross >> Max.P4 >> dross >> Steps[14] >> dross;
	is >> dross >> dross >> Min.Lambda >> dross >> Max.Lambda >> dross >> Steps[15] >> dross;
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
	vector<string> nv = (pr.Min).nameVector();
	vector<string> vvMin = (pr.Min).valueVector();
	vector<string> vvMax = (pr.Max).valueVector();
	for (uint j=0; j<pr.Size; j++)
		os << setw(20) << nv[j] << "[ " << setw(12) << vvMin[j] << " , " \
						<< setw(12) << vvMax[j] << " , " << setw(12) << (pr.Steps)[j] << " ]" << endl;
	return os;
}

// operator==
bool operator==(const ParametersRange& l, const ParametersRange& r) {
	return (l.Min==r.Min && l.Max==r.Max && l.Steps==r.Steps);
}

