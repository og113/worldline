/*
	definitions for program to generate unit loops in D dimensions for use in worldline programs.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h> 
#include "simple.h"
#include "genloop.h"
#include "genfloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - FCoeff class
		2 - functions acting on FCoeffs
		3 - FLoop class
		4 - FLoop functions
		5 - Metropolis
		6 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - FCoeff class
		- initializers and destructors
		- indexing
		- algebraic operations
----------------------------------------------------------------------------------------------------------------------------*/

// initializers and destructors
template<uint Dim>
FCoeff<Dim>::FCoeff(): Coeffs(2*Dim,0.0) {}

template <uint Dim>
FCoeff<Dim>::FCoeff(const FCoeff<Dim>& p): Coeffs(p.Coeffs) {}

template <uint Dim>
FCoeff<Dim>::FCoeff(const vector<number>& c): Coeffs(2*Dim) {
	if (c.size()==2*Dim)
		Coeffs = c;
	else
		cerr << "FCoeff error: cannot initialize Dim=" << Dim << " FCoeff from vector of size " << c.size() << endl;
}

template<uint Dim>
void FCoeff<Dim>::copy(const FCoeff<Dim>& p) {
	Coeffs = p.Coeffs;
}

template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator=(const FCoeff<Dim>& p) {
	copy(p);
	return *this;
}

template <uint Dim>
void FCoeff<Dim>::operator()(const FCoeff& p) {
	copy(p);
}

// set to zero
template<uint Dim>
void FCoeff<Dim>::zero() {
	fill(Coeffs.begin(), Coeffs.end(), 0.0);
}

// indexing
template <uint Dim>
const number& FCoeff<Dim>::operator[](const uint& loc) const {
	return Coeffs.at(loc);
}

template <uint Dim>
number& FCoeff<Dim>::operator[](const uint& loc) {
	return Coeffs.at(loc);
}

// +=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator+=(const FCoeff<Dim>& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] += rhs.Coeffs[j];
	return *this;
}

// -=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator-=(const FCoeff<Dim>& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] -= rhs.Coeffs[j];
	return *this;
}

// *=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator*=(const number& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] *= rhs;
	return *this;
}

// /=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator/=(const number& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] /= rhs;
	return *this;
}

// operator +
template <uint Dim>
FCoeff<Dim> operator+(const FCoeff<Dim>& lhs,const FCoeff<Dim>& rhs) {
	FCoeff<Dim> p(lhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] += rhs.Coeffs[j];
	return p;
}

// operator -
template <uint Dim>
FCoeff<Dim> operator-(const FCoeff<Dim>& lhs,const FCoeff<Dim>& rhs) {
	FCoeff<Dim> p(lhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] -= rhs.Coeffs[j];
	return p;
}

// operator *
template <uint Dim>
FCoeff<Dim> operator*(const number& lhs,const FCoeff<Dim>& rhs) {
	FCoeff<Dim> p(rhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] *= lhs;
	return p;
}

// operator /
template <uint Dim>
FCoeff<Dim> operator/(const FCoeff<Dim>& lhs,const number& rhs) {
	FCoeff<Dim> p(lhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] /= rhs;
	return p;
}


// operator<<
template <uint Dim>
ostream& operator<<(ostream& os,const FCoeff<Dim>& p) {
	os << left;
	for (uint j=0; j<2*Dim; j++)
		os << setw(15) << p[j];
	return os;
}

// operator>>
template <uint Dim>
istream& operator>>(istream& is,FCoeff<Dim>& p) {
	for (uint j=0; j<2*Dim; j++)
		is >> p[j];
	return is;
}

// writeBinary
template <uint Dim>
ostream& FCoeff<Dim>::writeBinary(ostream& os) const {
	for (uint j=0; j<2*Dim; j++)
		os.write(reinterpret_cast<const char*>(&(Coeffs[j])),sizeof(number));
	return os;
}

// readBinary
template <uint Dim>
istream& FCoeff<Dim>::readBinary(istream& is) {
	for (uint j=0; j<Dim; j++)
		is.read(reinterpret_cast<char*>(&(Coeffs[j])),sizeof(number));
	return is;
}	

// operator==
template <uint Dim>
bool operator==(const FCoeff<Dim>& lhs, const FCoeff<Dim>& rhs) {
	bool equal = true;
	for (uint j=0; j<2*Dim; j++) {
		if (abs(lhs[j]-rhs[j])>MIN_NUMBER) {
			equal = false;
			break;
		}
	}
	return equal;
}

// operator^=
template <uint Dim>
bool operator^=(const FCoeff<Dim>& lhs, const FCoeff<Dim>& rhs) {
	number closeness = MIN_NUMBER*1.0e4;
	bool equal = true;
	for (uint j=0; j<2*Dim; j++) {
		if (abs(lhs[j]-rhs[j])>closeness) {
			equal = false;
			break;
		}
	}
	return equal;
}


/*----------------------------------------------------------------------------------------------------------------------------
	2 - functions acting on FCoeffs

----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	3 - FLoop class
		- initialisers and destructor
		- steps
		- save
		- <<
----------------------------------------------------------------------------------------------------------------------------*/

// minimal initialiser
template <uint Dim>
Floop<Dim>::Floop(const uint& k, const uint& seed): 
		K(k), Seed(seed), Size(pow(2,k)), Grown(false) {
	FCoeffs.resize(Size);
	(FCoeffs[0]).zero();
	Generator = gsl_rng_alloc(gsl_rng_taus); // could also use gsl_rng_mt19937 (mersener twist)
}

// destructor
template <uint Dim>
Floop<Dim>::~Floop() {
	//gsl_rng_free(Generator);
	delete Generator;
}

// size
template <uint Dim>
uint Floop<Dim>::size() const {
	return Size;
}	

// clear
template <uint Dim>
void Floop<Dim>::clear() {
	for (uint j=0; j<Size; j++)
		FCoeffs[j].zero();
}

// set seed
template <uint Dim>
void Floop<Dim>::setSeed(const uint& s) {
	Seed = s;
}

// indexing
template <uint Dim>
const FCoeff<Dim>& Floop<Dim>::operator[](const uint& loc) const {
	return FCoeffs.at(loc);
}

// indexing
template <uint Dim>
FCoeff<Dim>& Floop<Dim>::operator[](const uint& loc) {
	return FCoeffs.at(loc);
}

// save
template <uint Dim>
void Floop<Dim>::save(const string& f) const {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (os.good()) {
		for (uint j=0; j<Size; j++) {
			(FCoeffs[j]).writeBinary(os);
		}
		os.close();
	}
	else {
		cerr << "Floop::save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// saveAscii
template <uint Dim>
void Floop<Dim>::saveAscii(const string& f) const {
	ofstream os;
	os.open(f.c_str());
	if (os.good()) {
		os << *this;
		os.close();
	}
	else {
		cerr << "Floop::save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// load
template <uint Dim>
void Floop<Dim>::load(const string& f) {
	if (!fileExists(f)) {
		cerr << "Floop::load error: file, " << f << ", doesn't exist" << endl;
		return;
	}
	uint fd = countDoubles(f);
	if (fd!=2*Dim*Size) {
		cerr << "Floop::load error: " << f << " contains " << fd << " doubles, floop requires " << 2*Size*Dim << endl;
		return;
	}
	else {	
		ifstream is;
		is.open(f.c_str(),ios::binary);
		if (is.good()) {
			for (uint j=0; j<Size; j++) {
				(FCoeffs[j]).readBinary(is);
			}
			is.close();
		}
		else {
			cerr << "Floop::load error: cannot read from " << f << endl;
			is.close();
			return;
		}
	}
	Grown = true;
}

// loadAscii
template <uint Dim>
void Floop<Dim>::loadAscii(const string& f) {
	uint fl = countLines(f);
	uint fc = countColumns(f);
	if (fl!=Size) {
		cerr << "Floop::load error: " << f << " contains " << fl << " lines, floop requires " << Size << endl;
		return;
	}
	else if (fc!=2*Dim) {
		cerr << "Floop::load error: " << f << " contains " << fc << " columns, floop requires " << 2*Dim << endl;
		return;
	}
	else {	
		ifstream is;
		is.open(f.c_str());
		if (is.good()) {
			is >> *this;
			is.close();
		}
		else {
			cerr << "Floop::load error: cannot read from " << f << endl;
			is.close();
			return;
		}
	}
	Grown = true;
}

// grow
template <uint Dim>
void Floop<Dim>::grow() {
	Generator = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(Generator,Seed);
	double sigma = SQRT2/pi;
		
	for (uint i=0; i<Size; i++) {
		for (uint j=0; j<Dim; j++) {
	  		(FCoeffs[i])[2*j] = gsl_ran_gaussian(Generator, sigma/(number)(i+1.0)); //gsl_ran_gaussian_ziggurat is another option
	  		(FCoeffs[i])[2*j+1] = gsl_ran_gaussian(Generator, sigma/(number)(i+1.0));
	  	}
	}

	//gsl_rng_free(Generator);
	Grown = true;
}

// getting spatial coord
template <uint Dim>
const Point<Dim> Floop<Dim>::X(const number& t) const {
	Point<Dim> x;
	double wt = 2.0*pi*t, wti;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		for (uint j=0; j<Dim; j++) {
			x[j] += (FCoeffs[i])[2*j]*cos(wti) + (FCoeffs[i])[2*j+1]*sin(wti);
		}
	}
	return x;
}

// getting spatial coord
template <uint Dim>
Point<Dim> Floop<Dim>::X(const number& t) {
	Point<Dim> x;
	double wt = 2.0*pi*t, wti;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		for (uint j=0; j<Dim; j++) {
			x[j] += (FCoeffs[i])[2*j]*cos(wti) + (FCoeffs[i])[2*j+1]*sin(wti);
		}
	}
	return x;
}


// stream <<
template <uint Dim>
ostream& operator<< (ostream& os,const Floop<Dim>& l) {
	for (uint j=0; j<l.Size; j++)
		os << l.FCoeffs[j] << endl;
	return os;
}

// stream >>
template <uint Dim>
istream& operator>> (istream& is, Floop<Dim>& l) {
	for (uint j=0; j<l.Size; j++)
		is >> l.FCoeffs[j];
	l.Grown = true;
	return is;
}

/*----------------------------------------------------------------------------------------------------------------------------
	4 - FLoop functions
		- S0
		- DS0
		- V0
		- aprxDV0
		
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	5 - FMetropolis
----------------------------------------------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------------------------------------------------------
	6 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------*/
// Dim=2
template class FCoeff<2>;
template ostream& operator<< <2>(ostream& os,const FCoeff<2>& p);
template FCoeff<2> operator+ <2>(const FCoeff<2>& lhs,const FCoeff<2>& rhs);
template FCoeff<2> operator- <2>(const FCoeff<2>& lhs,const FCoeff<2>& rhs);
template FCoeff<2> operator* <2>(const number& lhs,const FCoeff<2>& rhs);
template FCoeff<2> operator/ <2>(const FCoeff<2>& lhs,const number& rhs);
template bool operator== <2>(const FCoeff<2>& lhs, const FCoeff<2>& rhs);
template bool operator^= <2>(const FCoeff<2>& lhs, const FCoeff<2>& rhs);
template class Floop<2>;
template ostream& operator<< <2>(ostream& os,const Floop<2>& l);

// Dim=4
template class FCoeff<4>;
template ostream& operator<< <4>(ostream& os,const FCoeff<4>& p);
template FCoeff<4> operator+ <4>(const FCoeff<4>& lhs,const FCoeff<4>& rhs);
template FCoeff<4> operator- <4>(const FCoeff<4>& lhs,const FCoeff<4>& rhs);
template FCoeff<4> operator* <4>(const number& lhs,const FCoeff<4>& rhs);
template FCoeff<4> operator/ <4>(const FCoeff<4>& lhs,const number& rhs);
template bool operator== <4>(const FCoeff<4>& lhs, const FCoeff<4>& rhs);
template bool operator^= <4>(const FCoeff<4>& lhs, const FCoeff<4>& rhs);
template class Floop<4>;
template ostream& operator<< <4>(ostream& os,const Floop<4>& l);
