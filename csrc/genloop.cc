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

// operator>>
template <uint Dim>
istream& operator>>(istream& is,Point<Dim>& p) {
	for (uint j=0; j<Dim; j++)
		is >> p[j];
	return is;
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
	number closeness = MIN_NUMBER*1.0e4;
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
		- DistanceSquared
		- Distance
		- Dot (product)
----------------------------------------------------------------------------------------------------------------------------*/

//distance squared
template <uint Dim>
number DistanceSquared(const Point<Dim>& p1, const Point<Dim>& p2) {
	number d = pow(p1[0]-p2[0],2);
	for (uint j=1; j<Dim; j++) {
		d += pow(p1[j]-p2[j],2);
	}
	return d;
}


//distance
template <uint Dim>
number Distance(const Point<Dim>& p1, const Point<Dim>& p2) {
	return sqrt(DistanceSquared(p1,p2));
}

// Dot
template <uint Dim>
number Dot(const Point<Dim>& p1, const Point<Dim>& p2, const Point<Dim>& q1, const Point<Dim>& q2) {
	number d = (p1[0]-p2[0])*(q1[0]-q2[0]);
	for (uint j=1; j<Dim; j++) {
		d += (p1[j]-p2[j])*(q1[j]-q2[j]);
	}
	return d;
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
		K(k), Seed(seed), Length(pow(2,k)), Grown(false) {
	Points.resize(Length);
	(Points[0]).zero();
	Generator = gsl_rng_alloc(gsl_rng_taus); // could also use gsl_rng_mt19937 (mersener twist)
}

// destructor
template <uint Dim>
Loop<Dim>::~Loop() {
	//gsl_rng_free(Generator);
	delete[] Generator;
}

// size
template <uint Dim>
uint Loop<Dim>::size() const {
	return Length;
}	

// first step
template <uint Dim>
void Loop<Dim>::firstStep() {
	number sigma = 1.0/SQRT2;
	Points[Length/2] = Points[0];
	for (uint j=0; j<Dim; j++)
		(Points[Length/2])[j] += gsl_ran_gaussian (Generator, sigma); //gsl_ran_gaussian_ziggurat is another option
}

// following steps - problem when loc=N, instead of 0.
template <uint Dim>
void Loop<Dim>::followingSteps() {
	number sigma = 1.0/2.0;
	uint stepSize = 4;
	
	for (uint l=1; l<K; l++) {
		for (uint m=0; m<pow(2,l); m++) {
			uint locU = ((2*m+2)*Length/stepSize==Length? 0: (2*m+2)*Length/stepSize);
			Points[(2*m+1)*Length/stepSize] = (Points[locU]+Points[(2*m)*Length/stepSize])/2.0;
			for (uint n=0; n<Dim; n++) {
				(Points[(2*m+1)*Length/stepSize])[n] += gsl_ran_gaussian(Generator, sigma); //gsl_ran_gaussian_ziggurat is another option
			}
		}
		sigma /= SQRT2;
		stepSize *= 2;
	}
}

// metropolis
template <uint Dim>
void Loop<Dim>::metropolis() {
	if (!Grown) {
		cerr << "Loop error: cannot run metropolis before loop is grown" << endl;
		return;
	}
	// choosing location to change
	uint loc = (uint)(gsl_rng_uniform (Generator)*(size()+1.0));
	
	// calculating generation probability
	
	// generating new point according to generation probability
	
	// calculating acceptance probability
	
	// accepting new point according to acceptance probability
}

// normalise
template <uint Dim>
void Loop<Dim>::normalise() {
	number L = Distance(Points[Length-1],Points[0]);
	for (uint j=0; j<(Length-1); j++)
		L += Distance(Points[j+1],Points[j]);
	for (uint l=0; l<Length; l++)
		Points[l] /= L;
}

// centre
template <uint Dim>
void Loop<Dim>::centre() {
	Point<Dim> Xcm(Points[0]);
	for (uint j=1; j<Length; j++)
		Xcm += Points[j];
	for (uint l=0; l<Length; l++)
		Points[l] -= Xcm;
}

// clear
template <uint Dim>
void Loop<Dim>::clear() {
	for (uint j=0; j<Length; j++)
		Points[j].zero();
}

// set seed
template <uint Dim>
void Loop<Dim>::setSeed(const uint& s) {
	Seed = s;
}

// indexing
template <uint Dim>
const Point<Dim>& Loop<Dim>::operator[](const uint& loc) const {
	return Points.at(loc);
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

// saveAscii
template <uint Dim>
void Loop<Dim>::saveAscii(const string& f) const {
	ofstream os;
	os.open(f.c_str());
	if (os.good()) {
		os << *this;
		os.close();
	}
	else {
		cerr << "save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// load
template <uint Dim>
void Loop<Dim>::load(const string& f) {
	uint fd = countDoubles(f);
	if (fd!=Dim*Length) {
		cerr << "load error: " << f << " contains " << fd << " doubles, loop requires " << Length*Dim << endl;
		return;
	}
	else {	
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
	Grown = true;
}

// loadAscii
template <uint Dim>
void Loop<Dim>::loadAscii(const string& f) {
	uint fl = countLines(f);
	uint fc = countColumns(f);
	if (fl!=Length) {
		cerr << "load error: " << f << " contains " << fl << " lines, loop requires " << Length << endl;
		return;
	}
	else if (fc!=Dim) {
		cerr << "load error: " << f << " contains " << fc << " columns, loop requires " << Dim << endl;
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
			cerr << "load error: cannot read from " << f << endl;
			is.close();
			return;
		}
	}
	Grown = true;
}

// grow
template <uint Dim>
void Loop<Dim>::grow() {
	Generator = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(Generator,Seed);
	
	firstStep();
	followingSteps();
	//normalise();
	centre();
	
	//gsl_rng_free(Generator);
	Grown = true;
}

// checkLength
template <uint Dim>
number Loop<Dim>::checkLength() const {
	number L = Distance(Points[Length-1],Points[0]);
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

// stream >>
template <uint Dim>
istream& operator>> (istream& is, Loop<Dim>& l) {
	for (uint j=0; j<l.Length; j++)
		is >> l.Points[j];
	l.Grown = true;
	return is;
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

// Dim=2
template class Point<2>;
template ostream& operator<< <2>(ostream& os,const Point<2>& p);
template Point<2> operator+ <2>(const Point<2>& lhs,const Point<2>& rhs);
template Point<2> operator- <2>(const Point<2>& lhs,const Point<2>& rhs);
template Point<2> operator* <2>(const number& lhs,const Point<2>& rhs);
template Point<2> operator/ <2>(const Point<2>& lhs,const number& rhs);
template bool operator== <2>(const Point<2>& lhs, const Point<2>& rhs);
template bool operator^= <2>(const Point<2>& lhs, const Point<2>& rhs);
template number Distance(const Point<2>&, const Point<2>&);
template class Loop<2>;
template ostream& operator<< <2>(ostream& os,const Loop<2>& l);
