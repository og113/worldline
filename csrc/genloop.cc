/*
	definitions for program to generate unit loops in D dimensions for use in worldline programs.
*/

#include <iostream>
#include <iomanip>
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


// operator<<
template <uint Dim>
ostream& operator<<(ostream& os,const Point<Dim>& p) {
	os << left;
	for (uint j=0; j<Dim; j++)
		os << setw(15) << p[j];
	return os;
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
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	4 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------*/
// Dim=4
template class Point<4>;
template ostream& operator<< <4>(ostream& os,const Point<4>& p);
template Point<4> operator+ <4>(const Point<4>& lhs,const Point<4>& rhs);
template Point<4> operator- <4>(const Point<4>& lhs,const Point<4>& rhs);
template bool operator== <4>(const Point<4>& lhs, const Point<4>& rhs);
template bool operator^= <4>(const Point<4>& lhs, const Point<4>& rhs);
template number Distance(const Point<4>&, const Point<4>&);

