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
	4 - loop functions
		- S0
		- DS0
		- V0
		- aprxDV0
		
n.b. functions are defined for unit loops. for other loops the result must be multiplied by appropriate factors (depending on Dim) of T, the length of the loop.
----------------------------------------------------------------------------------------------------------------------------*/

// S0
template <uint Dim>
number S0 (const Loop<Dim>& l) {
	number result = DistanceSquared(l[l.size()-1],l[0]);
	for (uint j=0; j<(l.size()-1); j++)
		result += DistanceSquared(l[j+1],l[j]);
	return result*(number)l.size()/4.0;
}

// DS0
template <uint Dim>
number DS0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	uint loc_pos = (loc==(l.size()-1)? 0: loc);
	uint loc_neg = (loc==0? (l.size()-1): loc);
	Point<Dim> a = l[loc_neg]+l[loc_pos]-l[loc];
	number result = Dot(p,l[loc],p,a);
	return result*(number)l.size()/2.0;
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
	number result = 2.0*Dot(l[2],l[1],l[1],l[0])*pow(DistanceSquared(l[1],l[0]),(2.0-Dim)/2.0);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += 2.0*Dot(l[posj],l[j],l[posk],l[k])*pow(DistanceSquared(l[j],l[k]),(2.0-Dim)/2.0);
		}
	}
	return result*l.size()/pow(l.size()-1.0,2);
}

// V1r
template <uint Dim>
number V1r (const Loop<Dim>& l, const number& a) {
	number result = 2.0*Dot(l[2],l[1],l[1],l[0])*pow(DistanceSquared(l[1],l[0])+a*a,(2.0-Dim)/2.0);
	uint posj, posk;
	for (uint j=2; j<l.size(); j++) {
		posj = (j!=(l.size()-1)?j+1:0);
		for (uint k=0; k<j; k++) {
			posk = (k!=(l.size()-1)?k+1:0);
			result += 2.0*Dot(l[posj],l[j],l[posk],l[k])*pow(DistanceSquared(l[j],l[k])+a*a,(2.0-Dim)/2.0);
		}
	}
	return result*l.size()/pow(l.size()-1.0,2);
}

// DV1
template <uint Dim>
number DV1 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc) {
	number result = 0.0;
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	for (uint j=0; j<l.size(); j++) {
		if (j!=loc) {
			posj = (j!=(l.size()-1)?j+1:0);
			if (posj==loc)
				result += 2.0*Dot(p,l[j],l[posloc],p)*pow(DistanceSquared(l[j],p),(2.0-Dim)/2.0);
			else
				result += 2.0*Dot(l[posj],l[j],l[posloc],p)*pow(DistanceSquared(l[j],p),(2.0-Dim)/2.0);
			result -= 2.0*Dot(l[posj],l[j],l[posloc],l[loc])*pow(DistanceSquared(l[j],l[loc]),(2.0-Dim)/2.0);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// DV1r
template <uint Dim>
number DV1r (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc, const number& a) {
	number result = 0.0;
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);
	for (uint j=0; j<l.size(); j++) {
		if (j!=loc) {
			posj = (j!=(l.size()-1)?j+1:0);
			if (posj==loc)
				result += 2.0*Dot(p,l[j],l[posloc],p)*pow(DistanceSquared(l[j],p)+a*a,(2.0-Dim)/2.0);
			else
				result += 2.0*Dot(l[posj],l[j],l[posloc],p)*pow(DistanceSquared(l[j],p)+a*a,(2.0-Dim)/2.0);
			result -= 2.0*Dot(l[posj],l[j],l[posloc],l[loc])*pow(DistanceSquared(l[j],l[loc])+a*a,(2.0-Dim)/2.0);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// I0
template <uint Dim>
number I0 (const Loop<Dim>& l) {
	return 0;
}

/*----------------------------------------------------------------------------------------------------------------------------
	5 - Metropolis
----------------------------------------------------------------------------------------------------------------------------*/

// constructor
template <uint Dim>
Metropolis<Dim>::Metropolis(Loop<Dim>& L, const Parameters& p, const uint& s): Seed(s), Steps(0) {
	LoopPtr = &L;
	P = &p;
	G = P->G;
	SOld = 0.0;

	// setting generator
	Generator = gsl_rng_alloc(gsl_rng_taus);
}

// constructor
template <uint Dim>
Metropolis<Dim>::Metropolis(const Parameters& p, const uint& s): Seed(s), Steps(0) {
	LoopPtr = NULL;
	P = &p;
	G = P->G;
	SOld = 0.0;

	// setting generator
	Generator = gsl_rng_alloc(gsl_rng_taus);
}

// destructor
template <uint Dim>
Metropolis<Dim>::~Metropolis() {
	delete Generator;
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
uint Metropolis<Dim>::step(const uint& Num) {
	gsl_rng_set(Generator,Seed);
	// checking loop grown
	if (!LoopPtr->Grown) {
		cerr << "Metropolis error: cannot run Metropolis before loop is grown" << endl;
		return 0.0;
	}
	// counter
	number counter = 0.0;
	// initializing S0
	SOld = S0(*LoopPtr);
	if (abs(G)>MIN_NUMBER)
		SOld += G*V1r(*LoopPtr,(*P).Epsi);
		
	// defining some parameters
	number sigma = 1.0/sqrt((number)LoopPtr->size());
	uint loc, loc_pos, loc_neg;
	number acc_prob, rand;
	Point<Dim> temp;

	
	// starting metropolis loop
	for (uint j=0; j<Num; j++) {
	
		// choosing location to change
		loc = (uint)(gsl_rng_uniform (Generator)*LoopPtr->size());
	
		// calculating generation probability and generating new point
		loc_pos = (loc==(LoopPtr->size()-1)? 0: loc+1);
		loc_neg = (loc==0? LoopPtr->size()-1: loc-1);
		temp = (LoopPtr->Points[loc_pos] + LoopPtr->Points[loc_neg])/2.0;
		for (uint n=0; n<Dim; n++)
			temp[n] += gsl_ran_gaussian(Generator, sigma); //gsl_ran_gaussian_ziggurat is another option
	
		// calculating change in action
		SChange = DS0<Dim>(*LoopPtr, temp, loc);
		if (abs(G)>MIN_NUMBER)
			SChange += G*DV1r<Dim>(*LoopPtr, temp, loc,(*P).Epsi);	
	
		// accepting new point according to acceptance probability
		if (SChange<0.0) {
			LoopPtr->Points[loc] = temp;
			SOld += SChange;
			counter++;
		}
		else {
			acc_prob = gsl_sf_exp(-SChange);
			rand = gsl_rng_uniform(Generator);
			if (rand<acc_prob) {
				LoopPtr->Points[loc] = temp;
				SOld += SChange;
				counter++;
			}
		}
		
	}
	Steps += counter;
	LoopPtr->centre();
	return Steps;
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
template number Distance(const Point<4>&, const Point<4>&);
template class Loop<4>;
template ostream& operator<< <4>(ostream& os,const Loop<4>& l);
template number S0<4> (const Loop<4>& l);
template class Metropolis<4>;

// V, Dim=4, slightly changed for speed
template <> number V0 <4>(const Loop<4>& l) {
	number result = 2.0/DistanceSquared(l[1],l[0]);
	for (uint j=2; j<l.size(); j++) {
		for (uint k=0; k<j; k++) {
			result += 2.0/DistanceSquared(l[j],l[k]);
		}
	}
	return result/pow(l.size()-1.0,2);
}

// DV, Dim=4, slightly changed for speed
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

// I, for Dim=4
template <> number I0 <4>(const Loop<4>& l) {
uint posj = 1;
	number result = (l[0])[1]*((l[1])[2]-(l[0])[2]);
	for (uint j=1; j<l.size(); j++) {
		posj = (j==(l.size()-1)?0:j+1);
		result += (l[j])[1]*((l[posj])[2]-(l[j])[2]);
	}
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
	return result*l.size()/pow(l.size()-1.0,2);
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
	return result*l.size()/pow(l.size()-1.0,2);
}

// DV1
template <> number DV1<4> (const Loop<4>& l, const Point<4>& p, const uint& loc) {
	number result = 0.0;
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc) {
			posj = (j!=(l.size()-1)?j+1:0);
			if (posj==loc)
				result += 2.0*Dot(p,l[j],l[posloc],p)/(DistanceSquared(l[j],p));
			else
				result += 2.0*Dot(l[posj],l[j],l[posloc],p)/(DistanceSquared(l[j],p));
			result -= 2.0*Dot(l[posj],l[j],l[posloc],l[loc])/(DistanceSquared(l[j],l[loc]));
		}
	}
	return result/pow(l.size()-1.0,2);
}

// DV1r
template <> number DV1r<4> (const Loop<4>& l, const Point<4>& p, const uint& loc, const number& a) {
	number result = 0.0;
	uint posj, posloc = (loc!=(l.size()-1)?loc+1:0);

	for (uint j=0; j<l.size(); j++) {
		if (j!=loc) {
			posj = (j!=(l.size()-1)?j+1:0);
			if (posj==loc)
				result += 2.0*Dot(p,l[j],l[posloc],p)/(DistanceSquared(l[j],p)+a*a);
			else
				result += 2.0*Dot(l[posj],l[j],l[posloc],p)/(DistanceSquared(l[j],p)+a*a);
			result -= 2.0*Dot(l[posj],l[j],l[posloc],l[loc])/(DistanceSquared(l[j],l[loc])+a*a);
		}
	}
	return result/pow(l.size()-1.0,2);
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
template number Distance(const Point<2>&, const Point<2>&);
template class Loop<2>;
template ostream& operator<< <2>(ostream& os,const Loop<2>& l);
template number S0<2> (const Loop<2>& l);
template class Metropolis<2>;

// V0, Dim=2, logarithmic, GF(x,y) = log(|x-y|)/2/pi
template <> number V0 <2>(const Loop<2>& l) {
	number result = 2.0/DistanceSquared(l[1],l[0]);
	for (uint j=2; j<l.size(); j++) {
		for (uint k=0; k<j; k++) {
			result += gsl_sf_log(DistanceSquared(l[j],l[k]));
		}
	}
	return result/pow(l.size()-1.0,2);
}

// DV0, Dim=2, logarithmic, GF(x,y) = log(|x-y|)/2/pi
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

