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
		4 - Loop functions
		5 - Metropolis
		6 - explicit template instantiation
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

// posNeigh
uint posNeigh(const uint& j, const uint& N) {
	return (j==(N-1)? 0: j+1);
}

// negNeigh
uint negNeigh(const uint& j, const uint& N) {
	return (j==0? (N-1): j-1);
}

// posNeighDisjoint
uint posNeighDisjoint(const uint& j, const uint& N) {
	if (j==(N/2-1))
		return 0;
	else if (j==(N-1))
		return N/2;
	else
		return j+1;
}

// negNeighDisjoint
uint negNeighDisjoint(const uint& j, const uint& N) {
	if (j==0)
		return N/2-1;
	else if (j==(N/2))
		return N-1;
	else
		return j-1;
}

// oppNeigh
uint oppNeigh(const uint& j, const uint& N) {
	return N-1-j;
}

//distance squared
template <uint Dim>
number DistanceSquared(const Point<Dim>& p1, const Point<Dim>& p2) {
	number d = pow(p1[0]-p2[0],2);
	for (uint j=1; j<Dim; j++) {
		d += pow(p1[j]-p2[j],2);
	}
	return d;
}


// distance
template <uint Dim>
number Distance(const Point<Dim>& p1, const Point<Dim>& p2) {
	return sqrt(DistanceSquared(p1,p2));
}

// spatial distance squared
template <uint Dim>
number SpatialDistanceSquared(const Point<Dim>& p1, const Point<Dim>& p2) {
	number d = pow(p1[0]-p2[0],2);
	for (uint j=1; j<(Dim-1); j++) {
		d += pow(p1[j]-p2[j],2);
	}
	return d;
}

// spatial distance
template <uint Dim>
number SpatialDistance(const Point<Dim>& p1, const Point<Dim>& p2) {
	return sqrt(SpatialDistanceSquared(p1,p2));
}

// Disjoint Distance squared
template <uint Dim>
number DistanceSquaredDisjoint(const Point<Dim>& p1, const Point<Dim>& p2, const number& beta) {
	number d = pow(mod<number>(p1[Dim-1]-p2[Dim-1],-beta/2.0,beta/2.0),2);
	for (uint j=0; j<(Dim-1); j++) {
		d += pow(p1[j]-p2[j],2);
	}
	return d;
}

// Disjoint Distance
template <uint Dim>
number DistanceDisjoint(const Point<Dim>& p1, const Point<Dim>& p2, const number& beta) {
	return sqrt(DistanceSquaredDisjoint(p1,p2,beta));
}

// NormSquared
template <uint Dim>
number NormSquared(const Point<Dim>& p) {
	number d = pow(p[0],2);
	for (uint j=1; j<Dim; j++) {
		d += pow(p[j],2);
	}
	return d;
}

// Norm
template <uint Dim>
number Norm(const Point<Dim>& p) {
	return sqrt(NormSquared(p));
}

// Dot
template <uint Dim>
number Dot(const Point<Dim>& p1, const Point<Dim>& p2) {
	number d = p1[0]*p2[0];
	for (uint j=1; j<Dim; j++) {
		d += p1[j]*p2[j];
	}
	return d;
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

// DotDisjoint
template <uint Dim>
number DotDisjoint(const Point<Dim>& p1, const Point<Dim>& p2, const number& beta) {
	number d = mod<number>(p1[Dim-1],-beta/2.0,beta/2.0)*mod<number>(p2[Dim-1],-beta/2.0,beta/2.0);
	for (uint j=0; j<(Dim-1); j++) {
		d += p1[j]*p2[j];
	}
	return d;
}

// DotDisjoint
template <uint Dim>
number DotDisjoint(const Point<Dim>& p1, const Point<Dim>& p2, const Point<Dim>& q1, const Point<Dim>& q2, const number& beta) {
	number d = mod<number>(p1[Dim-1]-p2[Dim-1],-beta/2.0,beta/2.0)*mod<number>(q1[Dim-1]-q2[Dim-1],-beta/2.0,beta/2.0);
	for (uint j=0; j<(Dim-1); j++) {
		d += (p1[j]-p2[j])*(q1[j]-q2[j]);
	}
	return d;
}

// Angle
template <uint Dim>
number Angle(const Point<Dim>& p1, const Point<Dim>& p2, const Point<Dim>& p3) {
	return acos(Dot(p1,p2,p2,p3)/Distance(p1,p2)/Distance(p2,p3));
}

// MidpointDistance
template <uint Dim>
number MidpointDistanceSquared(const Point<Dim>& p1, const Point<Dim>& p2, const Point<Dim>& p3, const Point<Dim>& p4) {
	return 0.25*DistanceSquared(p1+p2,p3+p4);
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
	clear();
	Generator = gsl_rng_alloc(gsl_rng_taus); // could also use gsl_rng_mt19937 (mersener twist)
}

// destructor
template <uint Dim>
Loop<Dim>::~Loop() {
	//gsl_rng_free(Generator);
	delete Generator;
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
	Xcm /= (number)Length;
	for (uint l=0; l<Length; l++)
		Points[l] -= Xcm;
}

// clear
template <uint Dim>
void Loop<Dim>::clear() {
	for (uint j=0; j<Length; j++)
		(Points[j]).zero();
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

// indexing
template <uint Dim>
Point<Dim>& Loop<Dim>::operator[](const uint& loc) {
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
		for (uint j=0; j<size(); j++)
			os << Points[j] << endl;
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
	if (!fileExists(f)) {
		cerr << "Loop::load error: file, " << f << ", doesn't exist" << endl;
		return;
	}
	uint fd = countDoubles(f);
	if (fd<Dim*Length) {
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
			for (uint j=0; j<size(); j++)
				is >> Points[j];
			Grown = true;
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

// length
template <uint Dim>
number Loop<Dim>::length() const {
	number L = Distance(Points[Length-1],Points[0]);
	for (uint j=0; j<(Length-1); j++)
		L += Distance(Points[j+1],Points[j]);
	return L;
}

// setLength
template <uint Dim>
void Loop<Dim>::setLength(const number& L) {
	number oldLength = length();
	if (oldLength>MIN_NUMBER) {
		number ratio = L/oldLength;
		for (uint j=0; j<Length; j++)
			Points[j] *= ratio;
	}
	else {
		cerr << "Loop::setLength error: oldLength = 0.0 " << endl;
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	4 - loop functions
		- S0
		- DS0
		- V0
		- aprxDV0
		
n.b. functions are defined for unit loops. for other loops the result must be multiplied by appropriate factors (depending on Dim) of T, the length of the loop.
----------------------------------------------------------------------------------------------------------------------------*/

// Dot
template <uint Dim>
number Dot(const Loop<Dim>& loop, const uint& i, const uint& j, const uint& k, const uint& l) {
	return Dot(loop[i],loop[j],loop[k],loop[l]);
}

// Dot
template <uint Dim>
number Dot(const Loop<Dim>& loop, const uint& i, const uint& j) {
	uint ni = (i==0? (loop.size()-1): i-1);
	uint nj = (j==0? (loop.size()-1): j-1);
	return Dot(loop[i],loop[ni],loop[j],loop[nj]);
}

// MidpointDistance
template <uint Dim>
number MidpointDistanceSquared(const Loop<Dim>& loop, const uint& i, const uint& j) {
	uint pi = (i==loop.size()? 0: i+1);
	uint pj = (j==loop.size()? 0: j+1);
	return MidpointDistanceSquared(loop[pi],loop[i],loop[pj],loop[j]);
}

// DX
template <uint Dim>
number DX(const Loop<Dim>& loop, const uint& i, const uint& mu) {
	uint ni = negNeigh(i,loop.size());
	return (loop[i])[mu]-(loop[ni])[mu];
}

// DX
template <uint Dim>
number DX(const Loop<Dim>& loop, const uint& i, const uint& j, const uint& mu) {
	return (loop[i])[mu]-(loop[j])[mu];
}

// DXDisjoint
template <uint Dim>
number DXDisjoint(const Loop<Dim>& loop, const uint& i, const uint& mu, const number& beta) {
	uint ni = negNeighDisjoint(i,loop.size());
	if (mu==(Dim-1))
		return mod<number>((loop[i])[mu]-(loop[ni])[mu],-beta/2.0,beta/2.0);
	else
		return (loop[i])[mu]-(loop[ni])[mu];
}

// DXDisjoint
template <uint Dim>
number DXDisjoint(const Loop<Dim>& loop, const uint& i, const uint& j, const uint& mu, const number& beta) {
	if (mu==(Dim-1))
		return mod<number>((loop[i])[mu]-(loop[j])[mu],-beta/2.0,beta/2.0);
	else
		return (loop[i])[mu]-(loop[j])[mu];
}

// L
template <uint Dim>
number L (const Loop<Dim>& l) {
	number Len = Distance(l[l.size()-1],l[0]);
	for (uint j=0; j<(l.size()-1); j++)
		Len += Distance(l[j+1],l[j]);
	return Len;
}

// DL
template <uint Dim>
number DL (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	uint pos = (loc==(l.size()-1)? 0 : loc+1);
	uint neg = (loc==0? (l.size()-1) : loc-1);
	return (Distance(l[pos],p) + Distance(p,l[neg]) - Distance(l[pos],l[loc]) - Distance(l[loc],l[neg]));
}

// Sm
template <uint Dim>
number Sm (const Loop<Dim>& l) {
	uint pos, neg;
	Point<Dim> pp;
	number res = 0.0, temp;
	for (uint j=0; j<l.size(); j++) {
		pos = (j==(l.size()-1)? 0 : j+1);
		neg = (j==0? (l.size()-1) : j-1);
		pp = l[j];
		pp *= 2.0;
		pp -= l[neg];
		temp = Dot(l[pos],pp,l[pos],l[j]);
		temp /= DistanceSquared(l[pos],l[j]);
		res += temp;
	}
	return res/(number)l.size();
}

// Sm
template <uint Dim>
number Sm (const Loop<Dim>& l, const uint& m, const uint& n) {
	uint pos, neg;
	Point<Dim> pp;
	number res = 0.0, temp;
	for (uint j=0; j<l.size(); j++) {
		if (j!=m && j!=n) {
			pos = (j==(l.size()-1)? 0 : j+1);
			neg = (j==0? (l.size()-1) : j-1);
			pp = l[j];
			pp *= 2.0;
			pp -= l[neg];
			temp = Dot(l[pos],pp,l[pos],l[j]);
			temp /= DistanceSquared(l[pos],l[j]);
			res += temp;
		}
	}
	return res/((number)l.size()-2.0);
}

// DSm
template <uint Dim>
number DSm (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	uint pos1 = (loc==(l.size()-1)? 0 : loc+1);
	uint pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
	uint neg1 = (loc==0? (l.size()-1) : loc-1);
	uint neg2 = (neg1==0? (l.size()-1) : neg1-1);
	Point<Dim> pp = l[neg1];
	pp *= 2.0;
	pp -= l[neg2];
	number res = Dot(p,pp,p,l[neg1]);
	res /= DistanceSquared(p,l[neg1]);
	number temp = Dot(l[loc],pp,l[loc],l[neg1]);
	temp /= DistanceSquared(l[loc],l[neg1]);
	res -= temp;
	temp = Dot(p,l[loc],l[pos2],l[pos1]);
	temp /= DistanceSquared(l[pos2],l[pos1]);
	res += temp;
	pp = l[pos1];
	pp += l[neg1];
	pp /= 2.0;
	temp = Dot(pp,p,l[pos1],p);
	temp /= DistanceSquared(l[pos1],p);
	temp *= 2.0;
	res += temp;
	temp = Dot(pp,l[loc],l[pos1],l[loc]);
	temp /= DistanceSquared(l[pos1],l[loc]);
	temp *= 2.0;
	res -= temp;
	return res/(number)l.size();
}

// KG
template <uint Dim> // N.B. not sure this is correct. only checked for circle. different evaluation in nrloop.cc
number KG (const Loop<Dim>& l) {
	uint pos1, pos2;
	Point<Dim> pp;
	number res = 0.0, t_sqrd, t_dot_sqrd, t_t_dot;
	for (uint j=0; j<l.size(); j++) {
		pos1 = (j==(l.size()-1)? 0 : j+1);
		pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
		pp = l[pos1];
		pp *= 2.0;
		pp -= l[pos2];
		t_sqrd = DistanceSquared(l[pos1],l[j]);
		t_dot_sqrd = DistanceSquared(l[j],pp);
		t_t_dot = Dot(l[j],pp,l[pos1],l[j]);
		res += sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd);
	}
	return res; // multiply and divide by N for dt and for averaging respectively
}

// KGMax
template <uint Dim>
number KGMax (const Loop<Dim>& l) {
	uint pos1, pos2;
	Point<Dim> pp;
	number res = 0.0, t_sqrd, t_dot_sqrd, t_t_dot, temp;
	for (uint j=0; j<l.size(); j++) {
		pos1 = (j==(l.size()-1)? 0 : j+1);
		pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
		pp = l[pos1];
		pp *= 2.0;
		pp -= l[pos2];
		t_sqrd = DistanceSquared(l[pos1],l[j]);
		t_dot_sqrd = DistanceSquared(l[j],pp);
		t_t_dot = Dot(l[j],pp,l[pos1],l[j]);
		temp = sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
		res = (temp>res? temp: res);
	}
	return res;
}

// KGMax
template <uint Dim>
number KGMax (const Loop<Dim>& l, const uint& ex1, const uint& ex2) {
	uint pos1, pos2;
	Point<Dim> pp;
	number res = 0.0, t_sqrd, t_dot_sqrd, t_t_dot, temp;
	for (uint j=0; j<l.size(); j++) {
		if (j!=ex1 && j!=ex2) {
			pos1 = (j==(l.size()-1)? 0 : j+1);
			pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
			pp = l[pos1];
			pp *= 2.0;
			pp -= l[pos2];
			t_sqrd = DistanceSquared(l[pos1],l[j]);
			t_dot_sqrd = DistanceSquared(l[j],pp);
			t_t_dot = Dot(l[j],pp,l[pos1],l[j]);
			temp = sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
			res = (temp>res? temp: res);
		}
	}
	return res;
}

// S0
template <uint Dim>
number S0 (const Loop<Dim>& l) {
	number result = DistanceSquared(l[l.size()-1],l[0]);
	for (uint j=0; j<(l.size()-1); j++)
		result += DistanceSquared(l[j+1],l[j]);
	return result*(number)l.size()/4.0;
}

// S0Disjoint
template <uint Dim>
number S0Disjoint (const Loop<Dim>& l, const number& beta) {
	number result = 0.0;
	uint pj;
	for (uint j=0; j<l.size(); j++) {
		pj = posNeighDisjoint(j,l.size());
		result += DistanceSquaredDisjoint(l[pj],l[j],beta);
	}
	return result*(number)l.size()/4.0;
}

// DS0
template <uint Dim>
number DS0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	uint loc_pos = (loc==(l.size()-1)? 0: loc+1);
	uint loc_neg = (loc==0? (l.size()-1): loc-1);
	number result = DistanceSquared(p,l[loc_neg]);
	result += DistanceSquared(l[loc_pos],p);
	result -= DistanceSquared(l[loc],l[loc_neg]);
	result -= DistanceSquared(l[loc_pos],l[loc]);
	return result*(number)l.size()/4.0;
}

// V0
template <uint Dim>
number V0 (const Loop<Dim>& l) {
	number result = 2.0*pow(DistanceSquared(l[1],l[0]),(2.0-Dim)/2.0);
	for (uint j=2; j<l.size(); j++) {
		for (uint k=0; k<j; k++) {
			result += 2.0*pow(DistanceSquared(l[j],l[k]),(2.0-Dim)/2.0);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// DV0
template <uint Dim>
number DV0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	number result = 0.0;
	for (uint j=0; j<l.size(); j++) {
		if (j!=loc) {
			result += 2.0*pow(DistanceSquared(l[j],p),(2.0-Dim)/2.0);
			result -= 2.0*pow(DistanceSquared(l[j],l[loc]),(2.0-Dim)/2.0);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// aprxDV0
template <uint Dim>
number aprxDV0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	// COMPLETELY WRONG
	number result = 0.0;
	for (uint j=1; j<l.size(); j++) {
		for (uint k=0; k<j; k++) {
			result += 2.0*pow(DistanceSquared(l[j],l[k]),-Dim/2.0);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// V1
template <uint Dim>
number V1 (const Loop<Dim>& l) {
	number result = Dot(l[2],l[1],l[1],l[0])*pow(DistanceSquared(l[1],l[0]),(2.0-Dim)/2.0);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += Dot(l[posj],l[j],l[posk],l[k])*pow(DistanceSquared(l[j],l[k]),(2.0-Dim)/2.0);
		}
	}
	return 2.0*result;
}

// V1r
template <uint Dim>
number V1r (const Loop<Dim>& l, const number& a) {
	number result = Dot(l[2],l[1],l[1],l[0])*pow(DistanceSquared(l[1],l[0])+a*a,(2.0-Dim)/2.0);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += Dot(l[posj],l[j],l[posk],l[k])*pow(DistanceSquared(l[j],l[k])+a*a,(2.0-Dim)/2.0);
		}
	}
	return 2.0*result;
}

// V2r
template <uint Dim>
number V2r (const Loop<Dim>& l, const number& a) {
	number result = Dot(l[2],l[1],l[1],l[0])*pow(0.25*DistanceSquared(l[2]+l[1],l[1]+l[0])+a*a,(2.0-Dim)/2.0);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += Dot(l[posj],l[j],l[posk],l[k])*pow(0.25*DistanceSquared(l[posj]+l[j],l[posk]+l[k])+a*a,(2.0-Dim)/2.0);
		}
	}
	return 2.0*result;
}

// DV1
template <uint Dim>
number DV1 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	uint negloc = (loc!=0?loc-1:(l.size()-1));
	
	number result = Dot(l[posloc],p,p,l[negloc])*pow(DistanceSquared(l[negloc],p),(2.0-Dim)/2.0)\
						- Dot(l[posloc],l[loc],l[loc],l[negloc])*pow(DistanceSquared(l[negloc],l[loc]),(2.0-Dim)/2.0);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc && j!=negloc) {
			posj = (j!=(l.size()-1)?j+1:0);
			result += Dot(p,l[loc],l[posj],l[j])*pow(DistanceSquared(l[negloc],l[j]),(2.0-Dim)/2.0)\
					+ Dot(l[posloc],p,l[posj],l[j])*pow(DistanceSquared(p,l[j]),(2.0-Dim)/2.0)\
					- Dot(l[posloc],l[loc],l[posj],l[j])*pow(DistanceSquared(l[loc],l[j]),(2.0-Dim)/2.0);
		}
	}
	return 2.0*result;
}

// DV1r
template <uint Dim>
number DV1r (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc, const number& a) {
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	uint negloc = (loc!=0?loc-1:(l.size()-1));
	
	number result = Dot(l[posloc],p,p,l[negloc])*pow(DistanceSquared(l[negloc],p)+a*a,(2.0-Dim)/2.0)\
						- Dot(l[posloc],l[loc],l[loc],l[negloc])*pow(DistanceSquared(l[negloc],l[loc])+a*a,(2.0-Dim)/2.0);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc && j!=negloc) {
			posj = (j!=(l.size()-1)?j+1:0);
			result += Dot(p,l[loc],l[posj],l[j])*pow(DistanceSquared(l[negloc],l[j])+a*a,(2.0-Dim)/2.0)\
					+ Dot(l[posloc],p,l[posj],l[j])*pow(DistanceSquared(p,l[j])+a*a,(2.0-Dim)/2.0)\
					- Dot(l[posloc],l[loc],l[posj],l[j])*pow(DistanceSquared(l[loc],l[j])+a*a,(2.0-Dim)/2.0);
		}
	}
	return 2.0*result;
}

// DV2r
template <uint Dim>
number DV2r (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc, const number& a) {
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	uint negloc = (loc!=0?loc-1:(l.size()-1));
	
	number result = Dot(l[posloc],p,p,l[negloc])*pow(0.25*DistanceSquared(l[negloc]+p,p+l[posloc])+a*a,(2.0-Dim)/2.0)\
						- Dot(l[posloc],l[loc],l[loc],l[negloc])*pow(0.25*DistanceSquared(l[negloc]+l[loc],l[loc]+l[posloc])+a*a,(2.0-Dim)/2.0);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc && j!=negloc) {
			posj = (j!=(l.size()-1)?j+1:0);
			result += Dot(p,l[loc],l[posj],l[j])*pow(0.25*DistanceSquared(l[negloc]+p,l[j]+l[posj])+a*a,(2.0-Dim)/2.0)\
					- Dot(l[loc],l[negloc],l[posj],l[j])*pow(0.25*DistanceSquared(l[negloc]+l[loc],l[j]+l[posj])+a*a,(2.0-Dim)/2.0)\
					+ Dot(l[posloc],p,l[posj],l[j])*pow(0.25*DistanceSquared(p+l[posloc],l[j]+l[posj])+a*a,(2.0-Dim)/2.0)\
					- Dot(l[posloc],l[loc],l[posj],l[j])*pow(0.25*DistanceSquared(l[loc]+l[posloc],l[j]+l[posj])+a*a,(2.0-Dim)/2.0);
		}
	}
	return 2.0*result;
}

// I0
template <uint Dim>
number I0 (const Loop<Dim>& l) {
	return 0;
}

// DI0
template <uint Dim>
number DI0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	return 0;
}

// FGamma = gamma*cot(gamma)-1
template <uint Dim>
number FGamma (const Loop<Dim>& l) {
	number result = -(number)l.size(), cot_gamma, temp;
	uint posj = 1, negj = l.size()-1;
	for (uint j=0; j<l.size(); j++) {	
		temp = Dot(l[posj],l[j],l[j],l[negj]);
		cot_gamma = DistanceSquared(l[posj],l[j])*DistanceSquared(l[j],l[negj]) - temp*temp;
		cot_gamma = temp/sqrt(cot_gamma);
		result += cot_gamma*atan(1.0/cot_gamma);
		posj = ((j+1)!=(l.size()-1)?j+2:0);
		negj = j;
	}
	return result;
}

// DFGamma
template <uint Dim>
number DFGamma (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	uint j = (loc!=0?loc-1:(l.size()-1));
	uint nj = (j!=0?j-1:(l.size()-1));
	uint pj = loc;
	number temp, cot_gamma, result = 0.0;
	
	for (uint k=0; k<3; k++) {
		temp = Dot(l[pj],l[j],l[j],l[nj]);
		cot_gamma = DistanceSquared(l[pj],l[j])*DistanceSquared(l[j],l[nj]) - temp*temp;
		cot_gamma = temp/sqrt(cot_gamma);
		result -= cot_gamma*atan(1.0/cot_gamma);
		
		temp = Dot((pj==loc?p:l[pj]),(j==loc?p:l[j]),(j==loc?p:l[j]),(nj==loc?p:l[nj]));
		cot_gamma = DistanceSquared((pj==loc?p:l[pj]),(j==loc?p:l[j]))*DistanceSquared((j==loc?p:l[j]),(nj==loc?p:l[nj])) - temp*temp;
		cot_gamma = temp/sqrt(cot_gamma);
		result += cot_gamma*atan(1.0/cot_gamma);
		
		nj = (nj==(l.size()-1)? 0 : nj+1);
		j = (j==(l.size()-1)? 0 : j+1);
		pj = (pj==(l.size()-1)? 0 : pj+1);
	}
	return result;
}

// interpolate
template <uint Dim>
void interpolate (const Loop<Dim>& in, Loop<Dim>& out) {
	if (in.K==out.K) {
		out = in;
	}
	else if (in.K>out.K) {
		uint N_out = pow(2,out.K);
		uint N_diff = pow(2,in.K-out.K);
		for (uint j=0; j<N_out; j++) {
			out[j] = in[N_diff*j];
		}
	}
	else {
		uint N_in = pow(2,in.K);
		uint N_out = pow(2,out.K);
		uint N_diff = pow(2,out.K-in.K);
		uint n, posn, m;
		number b_m;
		for (uint j=0; j<N_out; j++) {
			n = j/N_diff;
			posn = (n!=(N_in-1)?n+1:0);
			m = j%N_diff;
			b_m = m/(number)N_diff;
			out[j] = (1.0-b_m)*in[n] + b_m*in[posn];
		}
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	5 - Metropolis
----------------------------------------------------------------------------------------------------------------------------*/

// constructor
template <uint Dim>
Metropolis<Dim>::Metropolis(Loop<Dim>& L, const Parameters& p, const uint& s): Seed(s), Steps(0), Data(), DataChange() {
	LoopPtr = &L;
	P = &p;

	// setting generator
	Generator = gsl_rng_alloc(gsl_rng_taus);
}

// constructor
template <uint Dim>
Metropolis<Dim>::Metropolis(const Parameters& p, const uint& s): Seed(s), Steps(0), Data(), DataChange() {
	LoopPtr = NULL;
	P = &p;

	// setting generator
	Generator = gsl_rng_alloc(gsl_rng_taus);
}

// destructor
template <uint Dim>
Metropolis<Dim>::~Metropolis() {
	delete Generator;
	//delete P;
	//delete LoopPtr;
}

// set loop
template <uint Dim>
void Metropolis<Dim>::setLoop(Loop<Dim>& L) {
	LoopPtr = &L;
}

// set seed
template <uint Dim>
void Metropolis<Dim>::setSeed(const uint& s) {
	Seed = s;
}

// Step
template <uint Dim>
uint Metropolis<Dim>::step(const uint& Num, const bool& firstStep, MetropolisData& data) {
	gsl_rng_set(Generator,Seed);
	// checking loop grown
	if (!LoopPtr->Grown) {
		cerr << "Metropolis error: cannot run Metropolis before loop is grown" << endl;
		return 0.0;
	}
	// counter
	number counter = 0.0, epsi = P->Epsi, g = P->G, b = P->B, t = P->T;
	uint size = LoopPtr->size();
	if (firstStep) {
		// initializing Data
		Data.L = L(*LoopPtr);
		Data.S0 = S0(*LoopPtr);
		Data.S = Data.S0;
		Data.Sm = Sm(*LoopPtr);
		if (abs(g)>MIN_NUMBER) {
			Data.V = V1r(*LoopPtr,epsi);
			Data.S += g*Data.V;
			Data.S -= (abs(epsi)>MIN_NUMBER? g*PI*Data.L/epsi: 0.0);
			Data.I0 = I0(*LoopPtr);
			Data.S -= g*b*t*Data.I0;
		}
	}
		
	// defining some parameters
	number sigma = 1.0/sqrt((number)size);
	uint loc, loc_pos, loc_neg;
	number acc_prob, rand;
	Point<Dim> temp;
	
	// starting metropolis loop
	for (uint j=0; j<Num; j++) {
	
		// choosing location to change
		loc = (uint)(gsl_rng_uniform (Generator)*size);
	
		// calculating generation probability and generating new point
		loc_pos = (loc==(size-1)? 0: loc+1);
		loc_neg = (loc==0? size-1: loc-1);
		temp = (LoopPtr->Points[loc_pos] + LoopPtr->Points[loc_neg])/2.0;
		for (uint n=0; n<Dim; n++)
			temp[n] += gsl_ran_gaussian(Generator, sigma); //gsl_ran_gaussian_ziggurat is another option
	
		// calculating change in action
		DataChange.S0 = DS0(*LoopPtr, temp, loc);
		DataChange.S = DataChange.S0;
		DataChange.Sm = DSm(*LoopPtr, temp, loc);
		if (abs(g)>MIN_NUMBER) {
			DataChange.V = DV1r<Dim>(*LoopPtr, temp, loc, epsi);	
			DataChange.S += g*DataChange.V;
			DataChange.L = DL(*LoopPtr, temp, loc);
			DataChange.S -= (abs(epsi)>MIN_NUMBER? g*PI*DataChange.L/epsi: 0.0);
			DataChange.I0 = DI0(*LoopPtr, temp, loc);
			DataChange.S -= g*b*t*DataChange.I0;
		}
	
		// accepting new point according to acceptance probability
		if (DataChange.S<0.0) {
			LoopPtr->Points[loc] = temp;
			Data.S += DataChange.S;
			Data.L += DataChange.L;
			Data.S0 += DataChange.S0;
			Data.V += DataChange.V;
			Data.I0 += DataChange.I0;
			Data.Sm += DataChange.Sm;
			counter++;
		}
		else if (DataChange.S<-LOG_MIN_NUMBER) {
			acc_prob = gsl_sf_exp(-DataChange.S);
			rand = gsl_rng_uniform(Generator);
			if (rand<acc_prob) {
				LoopPtr->Points[loc] = temp;
				Data.S += DataChange.S;
				Data.L += DataChange.L;
				Data.S0 += DataChange.S0;
				Data.V += DataChange.V;
				Data.I0 += DataChange.I0;
				Data.Sm += DataChange.Sm;
				counter++;
			}
		}
		
	}
	Steps += counter;
	LoopPtr->centre();
	data = Data;
	return Steps;
}


// Step
template <uint Dim>
uint Metropolis<Dim>::step(const uint& Num, const bool& firstStep) {
	MetropolisData dross;
	return step(Num,firstStep,dross);
}


/*----------------------------------------------------------------------------------------------------------------------------
	6 - explicit template instantiation
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
template number Norm<4>(const Point<4>&);
template number Distance<4>(const Point<4>&, const Point<4>&);
template number DistanceSquared<4>(const Point<4>&, const Point<4>&);
template number SpatialDistance<4>(const Point<4>&, const Point<4>&);
template number SpatialDistanceSquared<4>(const Point<4>&, const Point<4>&);
template number DistanceDisjoint<4>(const Point<4>&, const Point<4>&, const number&);
template number DistanceSquaredDisjoint<4>(const Point<4>&, const Point<4>&, const number&);
template number Dot<4>(const Point<4>&, const Point<4>&);
template number Dot<4>(const Point<4>&, const Point<4>&, const Point<4>&, const Point<4>&);
template number DotDisjoint<4>(const Point<4>&, const Point<4>&, const number&);
template number DotDisjoint<4>(const Point<4>&, const Point<4>&, const Point<4>&, const Point<4>&, const number&);
template number Angle<4>(const Point<4>&, const Point<4>&, const Point<4>&);
template number MidpointDistanceSquared<4>(const Point<4>&, const Point<4>&, const Point<4>&, const Point<4>&);
template number Dot<4>(const Loop<4>&, const uint&, const uint&, const uint&, const uint&);
template number Dot<4>(const Loop<4>&, const uint&, const uint&);
template number MidpointDistanceSquared<4>(const Loop<4>&, const uint&, const uint&);
template number DX<4>(const Loop<4>&, const uint&, const uint&);
template number DX<4>(const Loop<4>&, const uint&, const uint&, const uint&);
template number DXDisjoint<4>(const Loop<4>&, const uint&, const uint&, const number&);
template number DXDisjoint<4>(const Loop<4>&, const uint&, const uint&, const uint&, const number&);
template class Loop<4>;
template number L<4> (const Loop<4>& l);
template number DL<4> (const Loop<4>& l, const Point<4>& p, const uint& loc);
template number Sm<4> (const Loop<4>& l);
template number Sm<4> (const Loop<4>& l, const uint& m, const uint& n);
template number DSm<4> (const Loop<4>& l, const Point<4>& p, const uint& loc);
template number KG<4> (const Loop<4>& l);
template number KGMax<4> (const Loop<4>& l);
template number KGMax<4> (const Loop<4>& l, const uint& ex1, const uint& ex2);
template number S0<4> (const Loop<4>& l);
template number S0Disjoint<4> (const Loop<4>& l, const number&);
template number DS0<4> (const Loop<4>& l, const Point<4>& p, const uint& loc);
template number FGamma<4> (const Loop<4>& l);
template number DFGamma<4> (const Loop<4>& l, const Point<4>& p, const uint& loc);
template class Metropolis<4>;
template void interpolate<4>(const Loop<4>& in, Loop<4>& out);

// V0, Dim=4, slightly changed for speed
template <> number V0 <4>(const Loop<4>& l) {
	number result = 2.0/DistanceSquared(l[1],l[0]);
	for (uint j=2; j<l.size(); j++) {
		for (uint k=0; k<j; k++) {
			result += 2.0/DistanceSquared(l[j],l[k]);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// DV0, Dim=4, slightly changed for speed
template <> number DV0 <4>(const Loop<4>& l, const Point<4>& p, const uint& loc) {
	number result = 0.0;
	for (uint j=0; j<l.size(); j++) {
		if (j!=loc) {
			result += 2.0/DistanceSquared(l[j],p);
			result -= 2.0/DistanceSquared(l[j],l[loc]);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// aprxDV0, Dim=4, slightly changed for speed
template <> number aprxDV0 <4> (const Loop<4>& l, const Point<4>& p, const uint& loc) {
	// COMPLETELY WRONG
	number result = 0.0;
	for (uint j=1; j<l.size(); j++) {
		for (uint k=0; k<j; k++) {
			result += 2.0/pow(DistanceSquared(l[j],l[k]),2.0);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// I0, for Dim=4
template <> number I0 <4>(const Loop<4>& l) {
	uint posj = 1;
	number result = (l[0])[2]*((l[posj])[3]-(l[0])[3]);
	for (uint j=1; j<l.size(); j++) {
		posj = (j==(l.size()-1)?0:j+1);
		result += (l[j])[2]*((l[posj])[3]-(l[j])[3]);
	}
	return result;
}

// DI0, for Dim=4
template <> number DI0 <4>(const Loop<4>& l, const Point<4>& p, const uint& loc) {
	uint pos = (loc==(l.size()-1)? 0: loc+1);
	uint neg = (loc==0? (l.size()-1):loc-1);
	number result = p[2]*((l[pos])[3]-p[3]);
	result -= (l[loc])[2]*((l[pos])[3]-(l[loc])[3]);
	result += (l[neg])[2]*(p[3]-(l[loc])[3]);
	return result;
}

// V1
template <> number V1<4> (const Loop<4>& l) {
	number result = 2.0*Dot(l[2],l[1],l[1],l[0])/DistanceSquared(l[1],l[0]);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += 2.0*Dot(l[posj],l[j],l[posk],l[k])/DistanceSquared(l[j],l[k]);
		}
	}
	return result;
}

// V1r
template <> number V1r<4> (const Loop<4>& l, const number& a) {
	number result = 2.0*Dot(l[2],l[1],l[1],l[0])/(DistanceSquared(l[1],l[0])+a*a);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += 2.0*Dot(l[posj],l[j],l[posk],l[k])/(DistanceSquared(l[j],l[k])+a*a);
		}
	}
	return result;
}

// V2r
template <> number V2r<4> (const Loop<4>& l, const number& a) {
	number result = 2.0*Dot(l[2],l[1],l[1],l[0])/(0.25*DistanceSquared(l[2]+l[1],l[1]+l[0])+a*a);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += 2.0*Dot(l[posj],l[j],l[posk],l[k])/(0.25*DistanceSquared(l[posj]+l[j],l[posk]+l[k])+a*a);
		}
	}
	return result;
}

// DV1
template <> number DV1<4> (const Loop<4>& l, const Point<4>& p, const uint& loc) {
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	uint negloc = (loc!=0?loc-1:(l.size()-1));
	
	number result = Dot(l[posloc],p,p,l[negloc])/DistanceSquared(l[negloc],p)\
						- Dot(l[posloc],l[loc],l[loc],l[negloc])/DistanceSquared(l[negloc],l[loc]);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc && j!=negloc) {
			posj = (j!=(l.size()-1)?j+1:0);
			result += Dot(p,l[loc],l[posj],l[j])/DistanceSquared(l[negloc],l[j])\
					+ Dot(l[posloc],p,l[posj],l[j])/DistanceSquared(p,l[j])\
					- Dot(l[posloc],l[loc],l[posj],l[j])/DistanceSquared(l[loc],l[j]);
		}
	}
	return 2.0*result;
}

// DV1r
template <> number DV1r<4> (const Loop<4>& l, const Point<4>& p, const uint& loc, const number& a) {
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	uint negloc = (loc!=0?loc-1:(l.size()-1));
	number a2 = a*a;
	
	number result = Dot(l[posloc],p,p,l[negloc])/(DistanceSquared(l[negloc],p)+a2)\
						- Dot(l[posloc],l[loc],l[loc],l[negloc])/(DistanceSquared(l[negloc],l[loc])+a2);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc && j!=negloc) {
			posj = (j!=(l.size()-1)?j+1:0);
			result += Dot(p,l[loc],l[posj],l[j])/(DistanceSquared(l[negloc],l[j])+a2)\
					+ Dot(l[posloc],p,l[posj],l[j])/(DistanceSquared(p,l[j])+a2)\
					- Dot(l[posloc],l[loc],l[posj],l[j])/(DistanceSquared(l[loc],l[j])+a2);
		}
	}
	return 2.0*result;
}

// DV2r
template <> number DV2r<4> (const Loop<4>& l, const Point<4>& p, const uint& loc, const number& a) {
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	uint negloc = (loc!=0?loc-1:(l.size()-1));
	
	number result = Dot(l[posloc],p,p,l[negloc])/(0.25*DistanceSquared(l[negloc]+p,p+l[posloc])+a*a)\
						- Dot(l[posloc],l[loc],l[loc],l[negloc])/(0.25*DistanceSquared(l[negloc]+l[loc],l[loc]+l[posloc])+a*a);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc && j!=negloc) {
			posj = (j!=(l.size()-1)?j+1:0);
			result += Dot(p,l[negloc],l[posj],l[j])/(0.25*DistanceSquared(l[negloc]+p,l[j]+l[posj])+a*a)\
					- Dot(l[loc],l[negloc],l[posj],l[j])/(0.25*DistanceSquared(l[negloc]+l[loc],l[j]+l[posj])+a*a)\
					+ Dot(l[posloc],p,l[posj],l[j])/(0.25*DistanceSquared(p+l[posloc],l[j]+l[posj])+a*a)\
					- Dot(l[posloc],l[loc],l[posj],l[j])/(0.25*DistanceSquared(l[loc]+l[posloc],l[j]+l[posj])+a*a);
		}
	}
	return 2.0*result;
}

// Dim=2
template class Point<2>;
template ostream& operator<< <2>(ostream& os,const Point<2>& p);
template Point<2> operator+ <2>(const Point<2>& lhs,const Point<2>& rhs);
template Point<2> operator- <2>(const Point<2>& lhs,const Point<2>& rhs);
template Point<2> operator* <2>(const number& lhs,const Point<2>& rhs);
template Point<2> operator/ <2>(const Point<2>& lhs,const number& rhs);
template bool operator== <2>(const Point<2>& lhs, const Point<2>& rhs);
template bool operator^= <2>(const Point<2>& lhs, const Point<2>& rhs);
template number Distance<2>(const Point<2>&, const Point<2>&);
template number DistanceSquared<2>(const Point<2>&, const Point<2>&);
template number SpatialDistance<2>(const Point<2>&, const Point<2>&);
template number SpatialDistanceSquared<2>(const Point<2>&, const Point<2>&);
template number DistanceDisjoint<2>(const Point<2>&, const Point<2>&, const number&);
template number DistanceSquaredDisjoint<2>(const Point<2>&, const Point<2>&, const number&);
template number Norm<2>(const Point<2>&);
template number Dot<2>(const Point<2>&, const Point<2>&);
template number Dot<2>(const Point<2>&, const Point<2>&, const Point<2>&, const Point<2>&);
template number DotDisjoint<2>(const Point<2>&, const Point<2>&, const number&);
template number DotDisjoint<2>(const Point<2>&, const Point<2>&, const Point<2>&, const Point<2>&, const number&);
template number Angle<2>(const Point<2>&, const Point<2>&, const Point<2>&);
template number MidpointDistanceSquared<2>(const Point<2>&, const Point<2>&, const Point<2>&, const Point<2>&);
template number Dot<2>(const Loop<2>&, const uint&, const uint&, const uint&, const uint&);
template number Dot<2>(const Loop<2>&, const uint&, const uint&);
template number MidpointDistanceSquared<2>(const Loop<2>&, const uint&, const uint&);
template number DX<2>(const Loop<2>&, const uint&, const uint&);
template number DX<2>(const Loop<2>&, const uint&, const uint&, const uint&);
template number DXDisjoint<2>(const Loop<2>&, const uint&, const uint&, const number&);
template number DXDisjoint<2>(const Loop<2>&, const uint&, const uint&, const uint&, const number&);
template class Loop<2>;
template number L<2> (const Loop<2>& l);
template number DL<2> (const Loop<2>& l, const Point<2>& p, const uint& loc);
template number S0<2> (const Loop<2>& l);
template number S0Disjoint<2> (const Loop<2>& l, const number&);
template class Metropolis<2>;
template void interpolate<2>(const Loop<2>& in, Loop<2>& out);

// V0, Dim=2, logarithmic, GF(x,y) = log(|x-y|)/2/PI
template <> number V0 <2>(const Loop<2>& l) {
	number result = 2.0/DistanceSquared(l[1],l[0]);
	for (uint j=2; j<l.size(); j++) {
		for (uint k=0; k<j; k++) {
			result += gsl_sf_log(DistanceSquared(l[j],l[k]));
		}
	}
	return result/pow(l.size()-1.0,2);
}

// DV0, Dim=2, logarithmic, GF(x,y) = log(|x-y|)/2/PI
template <> number DV0 <2>(const Loop<2>& l, const Point<2>& p, const uint& loc) {
	number result = 0.0;
	for (uint j=0; j<l.size(); j++) {
		if (j!=loc) {
			result += gsl_sf_log(DistanceSquared(l[j],p));
			result -= gsl_sf_log(DistanceSquared(l[j],l[loc]));
		}
	}
	return result/pow(l.size()-1.0,2);
}

