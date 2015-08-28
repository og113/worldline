/*
	header file for to generate unit loops in D dimensions for use in worldline programs.
*/
 
#ifndef __GENLOOP_H_INCLUDED__
#define __GENLOOP_H_INCLUDED__

#include <vector>
#include <gsl/gsl_randist.h> 	// Distributions of random numbers
#include <gsl/gsl_sf_exp.h> 	// Exponential functions

using namespace std;

typedef double number;
typedef unsigned int uint;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - Point
		2 - functions acting on Points
		3 - Loop class
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	Point class
----------------------------------------------------------------------------------------------------------------------------*/
template <uint Dim>
class Point;

// operator+
template <uint Dim>
Point<Dim> operator+(const Point<Dim>&,const Point<Dim>&);

// operator-
template <uint Dim>
Point<Dim> operator-(const Point<Dim>&,const Point<Dim>&);

// operator<<
template <uint Dim>
ostream& operator<<(ostream&,const Point<Dim>&);

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
	
	// indexing
	const number& operator[](const uint&) const;
	number& operator[](const uint&);
	
	// algebraic operations
	Point& operator+=(const Point& rhs);
	Point& operator-=(const Point& rhs);
	friend Point operator+ <Dim>(const Point&,const Point&);
	friend Point operator- <Dim>(const Point&,const Point&);
	
	// boolean checks
	friend bool operator== <Dim>(const Point& lhs, const Point& rhs);
	friend bool operator^= <Dim>(const Point& lhs, const Point& rhs);
	
	// stream <<
	friend ostream& operator<< <Dim>(ostream&,const Point<Dim>&);
	
private:
	void copy(const Point&);
	vector<number>		Coords;
};

/*----------------------------------------------------------------------------------------------------------------------------
	functions acting on Points
----------------------------------------------------------------------------------------------------------------------------*/
template <uint Dim>
number Distance(const Point<Dim>&, const Point<Dim>&);


/*----------------------------------------------------------------------------------------------------------------------------
	Loop class
----------------------------------------------------------------------------------------------------------------------------*/

/*class Loop {
public:
	Loop(const uint& k,
private:
	uint			Dim;
	vector< Point<Dim> >	Points;
};*/

#endif // __GENLOOP_H_INCLUDED__
