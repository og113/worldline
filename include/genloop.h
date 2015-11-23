/*
	header file for to generate unit loops in D dimensions for use in worldline programs.
*/
 
#ifndef __GENLOOP_H_INCLUDED__
#define __GENLOOP_H_INCLUDED__

#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_randist.h> 	// Distributions of random numbers
#include <gsl/gsl_sf_exp.h> 	// Exponential functions
#include <gsl/gsl_sf_log.h>		// Logarithm functions
//#include <gsl/gsl_sf_gamma.h>   // Gamma functions
#include "simple.h"
#include "parameters.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - Point
		2 - functions acting on Points
		3 - Loop class
		4 - Loop functions
		5 - Metropolis
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - Point class
----------------------------------------------------------------------------------------------------------------------------*/
template <uint Dim>
class Point;

// operator+
template <uint Dim>
Point<Dim> operator+(const Point<Dim>&,const Point<Dim>&);

// operator-
template <uint Dim>
Point<Dim> operator-(const Point<Dim>&,const Point<Dim>&);

// operator*
template <uint Dim>
Point<Dim> operator*(const number&,const Point<Dim>&);

// operator/
template <uint Dim>
Point<Dim> operator/(const Point<Dim>&,const number&);

// operator<<
template <uint Dim>
ostream& operator<<(ostream&,const Point<Dim>&);

// operator>>
template <uint Dim>
istream& operator>>(istream&,Point<Dim>&);

// operator==
template <uint Dim>
bool operator==(const Point<Dim>& lhs, const Point<Dim>& rhs);

// operator^=
template <uint Dim>
bool operator^=(const Point<Dim>& lhs, const Point<Dim>& rhs);


// Point class
template <uint Dim>
class Point {
public:
	// initialization and destruction
	Point();
	Point(const Point&);
	Point(const vector<number>& c);
	~Point() {}
	
	// copying
	Point& operator=(const Point&);
	void operator()(const Point&);
	
	// set to zero
	void zero();
	
	// indexing
	const number& operator[](const uint&) const;
	number& operator[](const uint&);
	
	// algebraic operations
	Point& operator+=(const Point& rhs);
	Point& operator-=(const Point& rhs);
	Point& operator*=(const number& rhs);
	Point& operator/=(const number& rhs);
	friend Point operator+ <Dim>(const Point&,const Point&);
	friend Point operator- <Dim>(const Point&,const Point&);
	friend Point operator* <Dim>(const number&,const Point&);
	friend Point operator/ <Dim>(const Point&,const number&);
	
	// boolean checks
	friend bool operator== <Dim>(const Point& lhs, const Point& rhs);
	friend bool operator^= <Dim>(const Point& lhs, const Point& rhs);
	
	// stream <<
	friend ostream& operator<< <Dim>(ostream&,const Point&);
	friend istream& operator>> <Dim>(istream&,Point&);
	
	// read and write binary
	ostream& writeBinary(ostream& os) const;
	istream& readBinary(istream& is);	
	
private:
	void copy(const Point&);
	vector<number>		Coords;
};

/*----------------------------------------------------------------------------------------------------------------------------
	2 - functions acting on Points
----------------------------------------------------------------------------------------------------------------------------*/
// Distance squared
template <uint Dim>
number DistanceSquared(const Point<Dim>&, const Point<Dim>&);

// Distance
template <uint Dim>
number Distance(const Point<Dim>&, const Point<Dim>&);

// Dot
template <uint Dim>
number Dot(const Point<Dim>&, const Point<Dim>&);

// Dot
template <uint Dim>
number Dot(const Point<Dim>&, const Point<Dim>&, const Point<Dim>&, const Point<Dim>&);

/*----------------------------------------------------------------------------------------------------------------------------
	3 - Loop class
----------------------------------------------------------------------------------------------------------------------------*/
template <uint Dim>
class Loop;

// Metropolis
template <uint Dim>
class Metropolis;

// Loop class
template <uint Dim>
class Loop {
public:
	// initialization and destruction
	Loop(const uint& k, const uint& seed);
	~Loop();
	
	// grow loop
	void grow();
	
	// Metropolis
	friend class Metropolis<Dim>;
	
	// clear
	void clear();
	
	// set seed
	void setSeed(const uint&);
	
	// save
	void save(const string& file) const;
	void saveAscii(const string& file) const;
	
	// load
	void load(const string& file);
	void loadAscii(const string& file);
	
	// indexing
	const Point<Dim>& operator[](const uint&) const;
	Point<Dim>& operator[](const uint&);
	
	// size
	uint size() const;
	
	// check length
	number					length() const;
	void					setLength(const number&);
	
private:
	uint					K;
	uint					Seed;
	uint					Length;
	vector< Point<Dim> >	Points;
  	gsl_rng* 				Generator;
  	bool					Grown;
  	void					firstStep();
	void 					followingSteps();
	void					normalise();
	void 					centre();
};

/*----------------------------------------------------------------------------------------------------------------------------
	4 - loop functions
		- L
		- S0
		- DS0
		- V0
		- aprxDV0
----------------------------------------------------------------------------------------------------------------------------*/

// L
template <uint Dim>
number L (const Loop<Dim>& l);

// DL
template <uint Dim>
number DL (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// Sm
template <uint Dim>
number Sm (const Loop<Dim>& l);

// DSm
template <uint Dim>
number DSm (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// KG
template <uint Dim>
number KG (const Loop<Dim>& l);

// S0
template <uint Dim>
number S0 (const Loop<Dim>& l);

// DS0
template <uint Dim>
number DS0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// V0
template <uint Dim>
number V0 (const Loop<Dim>& l);

// DV0
template <uint Dim>
number DV0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// aprxDV0
template <uint Dim>
number aprxDV0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// V1
template <uint Dim>
number V1 (const Loop<Dim>& l);

// V1r
template <uint Dim>
number V1r (const Loop<Dim>& l, const number& epsi);

// DV1
template <uint Dim>
number DV1 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// DV1r
template <uint Dim>
number DV1r (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc, const number& a);

// I0
template <uint Dim>
number I0 (const Loop<Dim>& l);

// DI0
template <uint Dim>
number DI0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// FGamma = gamma*cot(gamma)-1
template <uint Dim>
number FGamma (const Loop<Dim>& l);

// DFGamma
template <uint Dim>
number DFGamma (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);


/*----------------------------------------------------------------------------------------------------------------------------
	5 - Metropolis
----------------------------------------------------------------------------------------------------------------------------*/

// Metropolis Data
struct MetropolisData {
	number L;
	number S0;
	number I0;
	number Vr;
	number FGamma;
	number S;
	number Sm;
};

// Metropolis
template <uint Dim>
class Metropolis {
public:
	Metropolis(Loop<Dim>& loop, const Parameters& p, const uint& seed);
	Metropolis(const Parameters& p, const uint& seed);
	~Metropolis();
	
	// step
	uint				step(const uint&, const bool& firstStep);
	uint				step(const uint&, const bool& firstStep, MetropolisData& data);
	
	// set seed, loop
	void 				setSeed(const uint&);
	void 				setLoop(Loop<Dim>&);
	
private:
	uint				Seed;
	uint				Steps;
	MetropolisData		Data;
	MetropolisData		DataChange;
	gsl_rng* 			Generator;
	const Parameters*		P;
	Loop<Dim>*			LoopPtr;
};

#endif // __GENLOOP_H_INCLUDED__
