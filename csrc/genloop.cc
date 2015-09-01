/*
	definitions for program to generate unit loops in D dimensions for use in worldline programs.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "simple.h"
#include "genloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - Point class
		2 - functions acting on Points
		3 - Loop class
		4 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - Point class
		- initializers and destructors
		- indexing
		- algebraic operations
----------------------------------------------------------------------------------------------------------------------------*/

// initializers and destructors
template<uint Dim>
Point<Dim>::Point(): Coords(Dim,0.0) {}

template <uint Dim>
Point<Dim>::Point(const Point<Dim>& p): Coords(p.Coords) {}

template <uint Dim>
Point<Dim>::Point(const vector<number>& c): Coords(Dim) {
	if (c.size()==Dim)
		Coords = c;
	else
		cerr << "Point error: cannot initialize Dim=" << Dim << " Point from vector of size " << c.size() << endl;
}

template<uint Dim>
void Point<Dim>::copy(const Point<Dim>& p) {
	Coords = p.Coords;
}

template <uint Dim>
Point<Dim>& Point<Dim>::operator=(const Point<Dim>& p) {
	copy(p);
	return *this;
}

template <uint Dim>
void Point<Dim>::operator()(const Point& p) {
	copy(p);
}

// set to zero
template<uint Dim>
void Point<Dim>::zero() {
	fill(Coords.begin(), Coords.end(), 0.0);
}

// indexing
template <uint Dim>
const number& Point<Dim>::operator[](const uint& loc) const {
	return Coords.at(loc);
}

template <uint Dim>
number& Point<Dim>::operator[](const uint& loc) {
	return Coords.at(loc);
}

// +=
template <uint Dim>
Point<Dim>& Point<Dim>::operator+=(const Point<Dim>& rhs) {
	for (uint j=0; j<Dim; j++)
		Coords[j] += rhs.Coords[j];
	return *this;
}

// -=
template <uint Dim>
Point<Dim>& Point<Dim>::operator-=(const Point<Dim>& rhs) {
	for (uint j=0; j<Dim; j++)
		Coords[j] -= rhs.Coords[j];
	return *this;
}

// *=
template <uint Dim>
Point<Dim>& Point<Dim>::operator*=(const number& rhs) {
	for (uint j=0; j<Dim; j++)
		Coords[j] *= rhs;
	return *this;
}

// /=
template <uint Dim>
Point<Dim>& Point<Dim>::operator/=(const number& rhs) {
	for (uint j=0; j<Dim; j++)
		Coords[j] /= rhs;
	return *this;
}

// operator +
template <uint Dim>
Point<Dim> operator+(const Point<Dim>& lhs,const Point<Dim>& rhs) {
	Point<Dim> p(lhs);
	for (uint j=0; j<Dim; j++)
		p.Coords[j] += rhs.Coords[j];
	return p;
}

// operator -
template <uint Dim>
Point<Dim> operator-(const Point<Dim>& lhs,const Point<Dim>& rhs) {
	Point<Dim> p(lhs);
	for (uint j=0; j<Dim; j++)
		p.Coords[j] -= rhs.Coords[j];
	return p;
}

// operator *
template <uint Dim>
Point<Dim> operator*(const number& lhs,const Point<Dim>& rhs) {
	Point<Dim> p(rhs);
	for (uint j=0; j<Dim; j++)
		p.Coords[j] *= lhs;
	return p;
}

// operator /
template <uint Dim>
Point<Dim> operator/(const Point<Dim>& lhs,const number& rhs) {
	Point<Dim> p(lhs);
	for (uint j=0; j<Dim; j++)
		p.Coords[j] /= rhs;
	return p;
}


// operator<<
template <uint Dim>
ostream& operator<<(ostream& os,const Point<Dim>& p) {
	os << left;
	for (uint j=0; j<Dim; j++)
		os << setw(15) << p[j];
	return os;
}

// writeBinary
template <uint Dim>
ostream& Point<Dim>::writeBinary(ostream& os) const {
	for (uint j=0; j<Dim; j++)
		os.write(reinterpret_cast<const char*>(&(Coords[j])),sizeof(number));
	return os;
}

// readBinary
template <uint Dim>
istream& Point<Dim>::readBinary(istream& is) {
	for (uint j=0; j<Dim; j++)
		is.read(reinterpret_cast<char*>(&(Coords[j])),sizeof(number));
	return is;
}	

// operator==
template <uint Dim>
bool operator==(const Point<Dim>& lhs, const Point<Dim>& rhs) {
	bool equal = true;
	for (uint j=0; j<Dim; j++) {
		if (abs(lhs[j]-rhs[j])>MIN_NUMBER) {
			equal = false;
			break;
		}
	}
	return equal;
}

// operator^=
template <uint Dim>
bool operator^=(const Point<Dim>& lhs, const Point<Dim>& rhs) {
	double closeness = MIN_NUMBER*1.0e4;
	bool equal = true;
	for (uint j=0; j<Dim; j++) {
		if (abs(lhs[j]-rhs[j])>closeness) {
			equal = false;
			break;
		}
	}
	return equal;
}


/*----------------------------------------------------------------------------------------------------------------------------
	2 - functions acting on Points
		- Distance
----------------------------------------------------------------------------------------------------------------------------*/

//distance
template <uint Dim>
number Distance(const Point<Dim>& p1, const Point<Dim>& p2) {
	number d = 0.0;
	for (uint j=0; j<Dim; j++) {
		d += pow(p1[j]-p2[j],2);
	}
	return sqrt(d);
}


/*----------------------------------------------------------------------------------------------------------------------------
	3 - Loop class
		- initialisers and destructor
		- steps
		- save
		- <<
----------------------------------------------------------------------------------------------------------------------------*/

// minimal initialiser
template <uint Dim>
Loop<Dim>::Loop(const uint& k, const uint& seed): 
		K(k), Seed(seed), Length(pow(2,k)) {
	Points.resize(Length);
	(Points[0]).zero();
	Generator = gsl_rng_alloc(gsl_rng_taus);
}

// destructor
template <uint Dim>
Loop<Dim>::~Loop() {
	//gsl_rng_free(Generator);
	delete[] Generator;
}	

// first step
template <uint Dim>
void Loop<Dim>::firstStep() {
	double sigma = 1.0/SQRT2;
	Points[Length/2] = Points[0];
	for (uint j=0; j<Dim; j++)
		(Points[Length/2])[j] += gsl_ran_gaussian_ziggurat (Generator, sigma);
}

// following steps - problem when loc=N, instead of 0.
template <uint Dim>
void Loop<Dim>::followingSteps() {
	double sigma = 1.0/2.0;
	uint stepSize = 4;
	
	for (uint l=1; l<K; l++) {
		for (uint m=0; m<pow(2,l); m++) {
			uint locU = ((2*m+2)*Length/stepSize==Length? 0: (2*m+2)*Length/stepSize);
			Points[(2*m+1)*Length/stepSize] = (Points[locU]+Points[(2*m)*Length/stepSize])/2.0;
			for (uint n=0; n<Dim; n++) {
				(Points[(2*m+1)*Length/stepSize])[n] += gsl_ran_gaussian_ziggurat(Generator, sigma);
			}
		}
		sigma /= SQRT2;
		stepSize *= 2;
	}
}

// normalise
template <uint Dim>
void Loop<Dim>::normalise() {
	double L = Distance(Points[Length-1],Points[0]);
	for (uint j=0; j<(Length-1); j++)
		L += Distance(Points[j+1],Points[j]);
	for (uint l=0; l<Length; l++)
		Points[l] /= L;
}

// clear
template <uint Dim>
void Loop<Dim>::clear() {
	for (uint j=0; j<(Length-1); j++)
		Points[j].zero();
}

// save
template <uint Dim>
void Loop<Dim>::save(const string& f) const {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (os.good()) {
		for (uint j=0; j<Length; j++) {
			(Points[j]).writeBinary(os);
		}
		os.close();
	}
	else {
		cerr << "save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// load - should add functionality to check length of file is correct
template <uint Dim>
void Loop<Dim>::load(const string& f) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	if (is.good()) {
		for (uint j=0; j<Length; j++) {
			(Points[j]).readBinary(is);
		}
		is.close();
	}
	else {
		cerr << "load error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

// grow
template <uint Dim>
void Loop<Dim>::grow() {
	Generator = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(Generator,Seed);
	
	firstStep();
	followingSteps();
	normalise();
	
	//gsl_rng_free(Generator);
}

// checkLength
template <uint Dim>
double Loop<Dim>::checkLength() const {
	double L = Distance(Points[Length-1],Points[0]);
	for (uint j=0; j<(Length-1); j++)
		L += Distance(Points[j+1],Points[j]);
	return L;
}

// stream <<
template <uint Dim>
ostream& operator<< (ostream& os,const Loop<Dim>& l) {
	for (uint j=0; j<l.Length; j++)
		os << l.Points[j] << endl;
	return os;
}

/*----------------------------------------------------------------------------------------------------------------------------
	4 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------*/
// Dim=4
template class Point<4>;
template ostream& operator<< <4>(ostream& os,const Point<4>& p);
template Point<4> operator+ <4>(const Point<4>& lhs,const Point<4>& rhs);
template Point<4> operator- <4>(const Point<4>& lhs,const Point<4>& rhs);
template Point<4> operator* <4>(const number& lhs,const Point<4>& rhs);
template Point<4> operator/ <4>(const Point<4>& lhs,const number& rhs);
template bool operator== <4>(const Point<4>& lhs, const Point<4>& rhs);
template bool operator^= <4>(const Point<4>& lhs, const Point<4>& rhs);
template number Distance(const Point<4>&, const Point<4>&);
template class Loop<4>;
template ostream& operator<< <4>(ostream& os,const Loop<4>& l);
