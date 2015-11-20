/*
	header file for to generate fourier representation of unit loops in D dimensions for use in worldline programs.
	
		should really have had these representations and the representation as spatial points being related by inheritance.
*/
 
#ifndef __GENFLOOP_H_INCLUDED__
#define __GENFLOOP_H_INCLUDED__

#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_randist.h> 	// Distributions of random numbers
#include <gsl/gsl_sf_exp.h> 	// Exponential functions
#include <gsl/gsl_sf_log.h>		// Logarithm functions
//#include <gsl/gsl_sf_gamma.h>   // Gamma functions
#include "simple.h"
#include "parameters.h"
#include "genloop.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - FCoeff
		2 - functions acting on FCoeffs
		3 - FLoop class
		4 - FLoop functions
		5 - FMetropolis
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - FCoeff class
----------------------------------------------------------------------------------------------------------------------------*/
template <uint Dim>
class FCoeff;

// operator+
template <uint Dim>
FCoeff<Dim> operator+(const FCoeff<Dim>&,const FCoeff<Dim>&);

// operator-
template <uint Dim>
FCoeff<Dim> operator-(const FCoeff<Dim>&,const FCoeff<Dim>&);

// operator*
template <uint Dim>
FCoeff<Dim> operator*(const number&,const FCoeff<Dim>&);

// operator/
template <uint Dim>
FCoeff<Dim> operator/(const FCoeff<Dim>&,const number&);

// operator<<
template <uint Dim>
ostream& operator<<(ostream&,const FCoeff<Dim>&);

// operator>>
template <uint Dim>
istream& operator>>(istream&,FCoeff<Dim>&);

// operator==
template <uint Dim>
bool operator==(const FCoeff<Dim>& lhs, const FCoeff<Dim>& rhs);

// operator^=
template <uint Dim>
bool operator^=(const FCoeff<Dim>& lhs, const FCoeff<Dim>& rhs);


// FCoeff class
template <uint Dim>
class FCoeff {
public:
	// initialization and destruction
	FCoeff();
	FCoeff(const FCoeff&);
	FCoeff(const vector<number>& a);
	~FCoeff() {}
	
	// copying
	FCoeff& operator=(const FCoeff&);
	void operator()(const FCoeff&);
	
	// set to zero
	void zero();
	
	// indexing
	const number& operator[](const uint&) const;
	number& operator[](const uint&);
	
	// algebraic operations
	FCoeff& operator+=(const FCoeff& rhs);
	FCoeff& operator-=(const FCoeff& rhs);
	FCoeff& operator*=(const number& rhs);
	FCoeff& operator/=(const number& rhs);
	friend FCoeff operator+ <Dim>(const FCoeff&,const FCoeff&);
	friend FCoeff operator- <Dim>(const FCoeff&,const FCoeff&);
	friend FCoeff operator* <Dim>(const number&,const FCoeff&);
	friend FCoeff operator/ <Dim>(const FCoeff&,const number&);
	
	// boolean checks
	friend bool operator== <Dim>(const FCoeff& lhs, const FCoeff& rhs);
	friend bool operator^= <Dim>(const FCoeff& lhs, const FCoeff& rhs);
	
	// stream <<
	friend ostream& operator<< <Dim>(ostream&,const FCoeff&);
	friend istream& operator>> <Dim>(istream&,FCoeff&);
	
	// read and write binary
	ostream& writeBinary(ostream& os) const;
	istream& readBinary(istream& is);	
	
private:
	void copy(const FCoeff&);
	vector<number>		Coeffs;
};

/*----------------------------------------------------------------------------------------------------------------------------
	2 - functions acting on FCoeffs
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	3 - FLoop class
----------------------------------------------------------------------------------------------------------------------------*/

template <uint Dim>
class Floop;

// stream <<
template <uint Dim>
ostream& operator<< (ostream&,const Floop<Dim>&);

// stream>>
template <uint Dim>
istream& operator>>(istream&,Floop<Dim>&);

// Metropolis
template <uint Dim>
class Metropolis;

// Floop class
template <uint Dim>
class Floop {
public:
	// initialization and destruction
	Floop(const uint& k, const uint& seed);
	~Floop();
	
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
	const FCoeff<Dim>& operator[](const uint&) const;
	FCoeff<Dim>& operator[](const uint&);
	
	// getting spatial coord
	const Point<Dim>& X(const number& t) const;
	Point<Dim>& X(const number& t);
	
	// stream <<, >>
	friend ostream& operator<< <Dim>(ostream&,const Floop&);
	friend istream& operator>> <Dim>(istream&, Floop&);
	
	// size
	uint size() const;
	
private:
	uint					K;
	uint					Seed;
	uint					Size;
	vector< FCoeff<Dim> >	FCoeffs;
  	gsl_rng* 				Generator;
  	bool					Grown;
};

/*----------------------------------------------------------------------------------------------------------------------------
	4 - FLoop functions	
----------------------------------------------------------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------------------------------------------------------
	5 - FMetropolis
----------------------------------------------------------------------------------------------------------------------------*/

// FMetropolis


#endif // __GENFLOOP_H_INCLUDED__
