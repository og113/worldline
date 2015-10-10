/*
	array
		declarations of some simple arrays
*/

#ifndef __CHECK_H_INCLUDED__
#define __CHECK_H_INCLUDED__

#include <iostream>
#include <string>
#include "simple.h"

/*-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. SimpleArray class
	2. MonteCarloData class
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. SimpleArray class
-------------------------------------------------------------------------------------------------------------------------*/

// declaration of SimpleArray
template <class T>
class SimpleArray;

// operator*
template <class T>
SimpleArray<T> operator*(const number&,const Point<Dim>&);

// operator/
template <class T>
SimpleArray<T> operator/(const number&,const Point<Dim>&);

// operator==
template <class T>
bool operator==(const SimpleArray<T>& lhs, const SimpleArray<T>& rhs);

// operator^=
template <class T>
bool operator^=(const SimpleArray<T>& lhs, const SimpleArray<T>& rhs);

template <class T>
class SimpleArray {
public:
	SimpleArray(const uint&);
	SimpleArray(const uint&, const T&);
	SimpleArray(const SimpleArray&);
	~SimpleArray();
	T& operator[](const uint&);
	const T& operator[](const uint&) const;
	
	// size
	uint size() const;
	
	// =, *=, /=
	SimpleArray& operator=(const SimpleArray&);
	SimpleArray& operator*=(const T&);
	SimpleArray& operator/=(const T&);
	
	// save and load
	void save(const string&) const;
	void load(const string&);

private:
	void copy(const SimpleArray&);
	uint Size;
	T* ArrayPtr;
};

/*-------------------------------------------------------------------------------------------------------------------------
	2. MonteCarloData class
-------------------------------------------------------------------------------------------------------------------------*/

class MonteCarloData {
public:
	// constructors, destructor
	MCData(const string&);
	~MCData();
	
	// save and load
	void save(const string&) const;
	void load(const string&);
	
	// calculations
	number 					calcMean();
	number 					calcMeanSqrd(); 
	SimpleArray<number>		calcCorrelator(); // don't really want to be copying something as large as a correlator around
	number 					calcIntCorrTime();
	number 					calcExpCorrTime();
	number 					calcVariance();
	number 					calcErrorJacknife();
	number 					calcErrorBootStrap(const uint&);
	
private:
	number Mean;
	number MeanSqrd;
	number IntCorrTime;
	number ExpCorrTime;
	number ErrorJacknife;
	number ErrorBootStrap;
	SimpleArray<number> DataArray;
	SimpleArray<number> Correlator;
};

#endif // __CHECK_H_INCLUDED__
