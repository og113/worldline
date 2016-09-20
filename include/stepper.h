/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for the Stepper class and dependencies
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __STEPPER_H_INCLUDED__
#define __STEPPER_H_INCLUDED__

#include <string>
#include <vector>
#include <utility> //for pair
#include <iostream>
#include "simple.h"
#include "error.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. Point2d
	2. Stepper
	3. static functions for testing
	
n.b. stepper defined in 2d
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. declarations for the simple Point2d class
		- Point2d
		- operator<<
		- operator+
		- operator-
		- operator^=
just a collection of two numbers (x,y)
-------------------------------------------------------------------------------------------------------------------------*/

// Point2d
class Point2d {
public:
	Point2d(): X(), Y() {}
	Point2d(const double& XX, const double& YY): X(XX), Y(YY) {}
	Point2d(const Point2d&);
	~Point2d() {}
	Point2d& operator=(const Point2d& rhs);
	void operator()(const double&, const double&);
	void operator()(const Point2d&);
	double X;
	double Y;
private:
	void copy(const Point2d&);
};

// operator<<
ostream& operator<<(ostream&,const Point2d&);

// operator>>
istream& operator>>(istream&,Point2d&);

// operator +
Point2d operator+(const Point2d&, const Point2d&);

// operator -
Point2d operator-(const Point2d&, const Point2d&);

// operator==
bool operator==(const Point2d& lhs, const Point2d& rhs);

// operator^=
bool operator^=(const Point2d& lhs, const Point2d& rhs);

/*-------------------------------------------------------------------------------------------------------------------------
	2. Stepper
		- FxyPair
		- StepperOptions
		- Stepper
-------------------------------------------------------------------------------------------------------------------------*/

// FxyPair
typedef pair<Point2d,double> FxyPair;

// <<
ostream& operator<<(ostream&, const FxyPair&);

// >>
istream& operator>>(istream&, FxyPair&);

// StepperOptions
struct StepperOptions{
	enum			stepTypeList {straight=1, constPlane=2};//, constTaylor=3, constLagrange=4};
	enum			directedList {undirected=1, global=2, local=3};
	double 			epsi_x;
	double			epsi_y;
	double 			angle0;
	double			tol;
	double			aim;
	bool			fixedAim;
	stepTypeList 	stepType;
	directedList	directed;
};

// <<
ostream& operator<<(ostream&, const StepperOptions&);

// >>
istream& operator>>(istream&, StepperOptions&);

// Stepper
class Stepper {
public:
	Stepper(const StepperOptions& sto, const double& X, const double& Y, const double& f0);
	Stepper(const StepperOptions& sto, const double& X, const double& Y);
	Stepper(const StepperOptions& sto, const Point2d& P, const double& f0);
	Stepper(const StepperOptions& sto, const Point2d& P);
	Stepper(const StepperOptions& sto);
	~Stepper() {}
	void 		setStart(const double& X, const double& Y);
	void 		setStart(const Point2d& P);
	void 		step();
	void 		addResult(const double& f);
	void		addResult(const double& f, const double& e, const double& n);
	void		save(const string&) const;
	void		load(const string&);
	Point2d 	point() const;
	Point2d 	point(const uint&) const;
	Point2d 	lastStep() const;
	double		x() const;
	double		y() const;
	double		result() const;
	double		result(const uint&) const;
	double		aim() const;
	uint		local() const;
	uint		steps() const;
	bool		keep() const;
	double		stepAngle() const;
	double		tol() const;
	
	//friend ostream& operator<<(ostream&, const Stepper&);
private:
	StepperOptions 	opts;
	vector<FxyPair> f_xy_local;
	vector<FxyPair> f_xy_steps;
	double			angle;
};

// <<
ostream& operator<<(ostream&, const Stepper&);

/*-------------------------------------------------------------------------------------------------------------------------
	3. static functions for testing
		- find_nth_closest
	n.b. to be removed once tested
-------------------------------------------------------------------------------------------------------------------------*/

// find_nth_closest
//uint find_nth_closest(const vector<FxyPair>& fxy, const double& f, const uint& n);

#endif // __STEPPER_H_INCLUDED__
