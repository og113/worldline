/*
	definitions for solving classical equations of motion for monopole worldline via newton-raphson method
*/

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <jbnumlib.h> // for the jacobi theta functions
#include "simple.h"
#include "parameters.h"
#include "genloop.h"
#include "nrloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - some static functions
		1 - nr loop functions
		2 - loopToVector, vectorToLoop
		3 - filename functions
		4 - explicit instatiation, for specific dimensions
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0 - some static functions
----------------------------------------------------------------------------------------------------------------------------*/

// FThermal
static number FThermal(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return  sinh((2.0*PI*ra)/beta)/(4.0*PI*ra*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))); //
}

// DFThermalDr - use Onr version instead as finite at r=0
/*static number DFThermalDr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (r*cosh((2.0*PI*ra)/beta))/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 - (r*sinh((2.0*PI*ra)/beta))/\
 (4.0*PI*pow(ra,2*1.5)*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + (r*pow(sinh((2.0*PI*ra)/beta),2))/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2));
}*/

// DFThermalDrOnr
static number DFThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return cosh((2.0*PI*ra)/beta)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 - sinh((2.0*PI*ra)/beta)/\
 (4.0*PI*pow(ra,3)*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + pow(sinh((2.0*PI*ra)/beta),2)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2));
 }

// DFThermalDt
static number DFThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (sin((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (2.0*ra*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2));
}
// DDFThermalDtDr - use Onr version instead as finite at r=0
/*static number DDFThermalDtDr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (PI*r*cosh((2.0*PI*ra)/beta)*sin((2.0*PI*t)/beta))/\
 ((pow(a,2) + pow(r,2))*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 - (r*sin((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (2.0*pow(ra,3)*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2))\
 + (2.0*PI*r*sin((2.0*PI*t)/beta)*pow(sinh((2.0*PI*ra)/beta),2))/\
 ((pow(a,2) + pow(r,2))*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),3));
}*/

// DDFThermalDrDr
static number DDFThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (-3.0*pow(r,2)*cosh((2.0*PI*ra)/beta))/\
 (2.0*pow(ra,4)*pow(beta,2)*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + cosh((2.0*PI*ra)/beta)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + (PI*pow(r,2)*sinh((2.0*PI*ra)/beta))/\
 (pow(ra,3)*pow(beta,3)*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + (3.0*pow(r,2)*sinh((2.0*PI*ra)/beta))/\
 (4.0*PI*pow(ra,5)*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 - sinh((2.0*PI*ra)/beta)/\
 (4.0*PI*pow(ra,3)*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + (3.0*PI*pow(r,2)*cosh((2.0*PI*ra)/beta)*sinh((2.0*PI*ra)/beta))/\
 (pow(ra,3)*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 - (3.0*pow(r,2)*pow(sinh((2.0*PI*ra)/beta),2))/\
 (2.0*pow(ra,4)*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 + pow(sinh((2.0*PI*ra)/beta),2)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 + (2.0*PI*pow(r,2)*pow(sinh((2.0*PI*ra)/beta),3))/\
 (pow(ra,3)*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),3));
}

// DDFThermalDtDr - use Onr version instead as finite at r=0
/*static number DDFThermalDtDr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*pow(PI,2)*sin((2.0*PI*t)/beta)*((2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta))*
 coshsinhc((2.0*PI*r)/beta) - 4.0*PI*r*beta*pow(sinch((2.0*PI*r)/beta),2)))/\
 pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}*/

 // DDFThermalDtDrOnr
static number DDFThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (PI*cosh((2.0*PI*ra)/beta)*sin((2.0*PI*t)/beta))/\
 ((pow(a,2) + pow(r,2))*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 - (sin((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (2.0*pow(ra,2*1.5)*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2))\
 + (2.0*PI*sin((2.0*PI*t)/beta)*pow(sinh((2.0*PI*ra)/beta),2))/\
 ((pow(a,2) + pow(r,2))*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),3));
}

// DDFThermalDtDt
static number DDFThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return  (PI*cos((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (ra*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 + (2.0*PI*pow(sin((2.0*PI*t)/beta),2)*sinh((2.0*PI*ra)/beta))/\
 (ra*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),3));
}

// GThermal
static number GThermal(const number& r, const number& t, const number& beta, const number& a) {
	return (a*sqrt(PI)*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/(exp(pow(r,2)/pow(a,2))*beta);
}

// DGThermalDr - use Onr version instead as finite at r=0
/*static number DGThermalDr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*r*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/(a*exp(pow(r,2)/pow(a,2))*beta);
}*/

// DGThermalDrOnr
static number DGThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/(a*exp(pow(r,2)/pow(a,2))*beta);
}

// DGThermalDt
static number DGThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	return (a*sqrt(PI)*jbderivutheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (exp(pow(r,2)/pow(a,2))*pow(beta,2));
}

// DDGThermalDrDr
static number DDGThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*(pow(a,2) - 2.0*pow(r,2))*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (pow(a,3)*exp(pow(r,2)/pow(a,2))*beta);
}

// DDGThermalDtDr - use Onr version instead as finite at r=0
/*static number DDGThermalDtDr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*r*jbderivutheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (a*exp(pow(r,2)/pow(a,2))*pow(beta,2));
}*/

// DDGThermalDtDrOnr
static number DDGThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*jbderivutheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (a*exp(pow(r,2)/pow(a,2))*pow(beta,2));
}

// DDGThermalDtDt
static number DDGThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	return (a*sqrt(PI)*jbderiv2utheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (exp(pow(r,2)/pow(a,2))*pow(beta,3));
}

// HThermal
static number HThermal(const number& r, const number& t, const number& beta, const number& a) {
	return  exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)));
}

// DHThermalDrOnr
static number DHThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2))))/pow(a,2);
}

// DHThermalDt
static number DHThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	return - ((exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*beta*sin((2.0*PI*t)/beta))/\
 (pow(a,2)*PI));
}

// DDHThermalDrDr
static number DDHThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2))))/pow(a,2) \
 + (4.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*pow(r,2))/pow(a,4);
}

// DDHThermalDtDrOnr
static number DDHThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*beta*sin((2.0*PI*t)/beta))/\
 (pow(a,4)*PI);
}

// DDHThermalDtDt
static number DDHThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	return  (-2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*cos((2.0*PI*t)/beta))/\
 pow(a,2) + (exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*pow(beta,2)*\
 pow(sin((2.0*PI*t)/beta),2))/(pow(a,4)*pow(PI,2));
}

/*----------------------------------------------------------------------------------------------------------------------------
	1 - nr loop functions
----------------------------------------------------------------------------------------------------------------------------*/

// L
template <uint Dim>
void L (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pos = posNeigh(j,l.size());
	result += f*Distance(l[pos],l[j]);
}

// LDisjoint
template <uint Dim>
void LDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pos = posNeighDisjoint(j,l.size());
	result += f*DistanceDisjoint(l[pos],l[j],beta);
}

// DistPow
template <uint Dim>
void DistPow (const uint& j, const Loop<Dim>& l, const number& w, const number& f, number& result) {
	uint pos = posNeigh(j,l.size());
	result += f*pow(DistanceSquared(l[pos],l[j]),w/2.0);
}

// Sm
template <uint Dim>
void Sm (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	Point<Dim> pp;
	pp = l[j];
	pp *= 2.0;
	pp -= l[nj];
	number temp = Dot(l[pj],pp,l[pj],l[j]);
	temp /= DistanceSquared(l[pj],l[j]);
	result += f*temp/(number)l.size();
}

// S0
template <uint Dim>
void S0 (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	result += f*DistanceSquared(l[pj],l[j])*(number)l.size()/4.0;
}


// S0Disjoint
template <uint Dim>
void S0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	result += f*DistanceSquaredDisjoint(l[pj],l[j],beta)*(number)l.size()/4.0;
}

// I0
template <uint Dim>
void I0 (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	result += f*(l[j])[Dim-2]*DX(l,pj,j,Dim-1);
}

// I0Disjoint
template <uint Dim>
void I0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	result += f*(l[j])[Dim-2]*DXDisjoint(l,pj,j,Dim-1,beta);
}

// In
template <uint Dim>
void In (const uint& j, const Loop<Dim>& l, const number& n, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	result += f*pow((l[j])[Dim-2],n)*DX(l,pj,j,Dim-1);
}

// InDisjoint
template <uint Dim>
void InDisjoint (const uint& j, const Loop<Dim>& l, const number& n, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	result += f*pow((l[j])[Dim-2],n)*DXDisjoint(l,pj,j,Dim-1,beta);
}

// Angle
template <uint Dim>
void Angle (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	number cos_gamma;
	
	cos_gamma = Dot(l[pj],l[j],l[j],l[nj]);
	cos_gamma /= Distance(l[pj],l[j])*Distance(l[j],l[nj]);
	
	result += f*(acos(cos_gamma));
}

// PseudoAngle, in 3-4 plane
template <uint Dim>
void PseudoAngle (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	cerr << "PseudoAngle error: not defined in dimension " << Dim << endl;
}

// FGamma
template <uint Dim>
void FGamma (const uint& j,const Loop<Dim>& l,  const number& f, number& result) {

	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	number cot_gamma, temp;
	
	temp = Dot(l[pj],l[j],l[j],l[nj]);
	cot_gamma = DistanceSquared(l[pj],l[j])*DistanceSquared(l[j],l[nj]) - temp*temp;
	cot_gamma = temp/sqrt(cot_gamma);
	
	result += f*(cot_gamma*atan(1.0/cot_gamma)-1.0);	
}

// MaxXn
template <uint Dim>
void MaxXn (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result) {
	if (n>=Dim) {
		cerr << "MaxXn error: " << n << ">=" << Dim << endl;
	}
	if (k<j) {
		number temp = abs((l[j])[n]-(l[k])[n]);
		result = (temp*f>result? temp*f: result);
	}
}

// MinXnDisjoint
template <uint Dim>
void MinXnDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result) {
	if (n>=Dim) {
		cerr << "MinXnDisjoint error: " << n << ">=" << Dim << endl;
	}
	if (k<(l.size()/2) && j>=(l.size()/2)) {
		number temp = abs((l[j])[n]-(l[k])[n]);
		result = (temp*f<result? temp*f: result);
	}
}

// Vor
template <uint Dim>
void Vor (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*pow(DistanceSquared(l[j],l[k])+a*a,(2.0-Dim)/2.0);	
	}
	
}

// Vlr
template <uint Dim>
void Vlr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*pow(0.25*DistanceSquared(l[pj]+l[j],l[pk]+l[k])+a*a,(2.0-Dim)/2.0);	
	}
}

// Ver
template <uint Dim>
void Ver (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
		number d2 = DistanceSquared(l[j],l[k]);
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])*pow(d2,(2.0-Dim)/2.0)*exp(-a*a/d2);
	}
}

// Vdr
template <uint Dim>
void Vdr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = posNeigh(j,l.size());
		uint pk = posNeigh(k,l.size());
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])*pow(DistanceSquared(l[j],l[k]),(2.0-Dim+a)/2.0);	
	}
}

// Vthr
template <uint Dim>
void Vthr (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	cerr << "Vthr error: not defined in dimension " << Dim << endl;
}

// VthrDisjoint
template <uint Dim>
void VthrDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	cerr << "VthrDisjoint error: not defined in dimension " << Dim << endl;
}

// VthrDisjointLR
template <uint Dim>
void VthrDisjointLR (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	cerr << "VthrDisjoint errorLR: not defined in dimension " << Dim << endl;
}

// VnonrelDisjoint
template <uint Dim>
void VnonrelDisjoint (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	cerr << "VnonrelDisjoint error: not defined in dimension " << Dim << endl;
}

// VnonrelrDisjoint
template <uint Dim>
void VnonrelrDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	cerr << "VnonrelrDisjoint error: not defined in dimension " << Dim << endl;
}

// Gaussian
template <uint Dim>
void Gaussian (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeigh(j,l.size());
		uint pk = posNeigh(k,l.size());
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*exp(-DistanceSquared(l[j],l[k])/a/a);	
	}
}

// GaussianDisjoint
template <uint Dim>
void GaussianDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeighDisjoint(j,l.size());
		uint pk = posNeighDisjoint(k,l.size());
	
		result += f*(1.0+(number)(k<j))*DotDisjoint(l[pj],l[j],l[pk],l[k],beta)*exp(-DistanceSquaredDisjoint(l[j],l[k],beta)/a/a);	
	}
}

// GaussianLRDisjoint
template <uint Dim>
void GaussianLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		GaussianDisjoint(j,k,l,beta,a,f,result);
}

// GaussianThermal
template <uint Dim>
void GaussianThermal (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeigh(j,l.size());
		uint pk = posNeigh(k,l.size());
		number r = SpatialDistance(l[j],l[k]);
		number t = DX(l,k,j,3);
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*GThermal(r,t,beta,a);	
	}
}

// GaussianThermal2
template <uint Dim>
void GaussianThermal2 (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeigh(j,l.size());
		uint pk = posNeigh(k,l.size());
		number r = SpatialDistance(l[j],l[k]);
		number t = DX(l,k,j,3);
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*HThermal(r,t,beta,a);	
	}
}

// GaussianThermalDisjoint
template <uint Dim>
void GaussianThermalDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeighDisjoint(j,l.size());
		uint pk = posNeighDisjoint(k,l.size());
		number r = SpatialDistance(l[j],l[k]);
		number t = DXDisjoint(l,k,j,3,beta);
	
		result += f*(1.0+(number)(k<j))*DotDisjoint(l[pj],l[j],l[pk],l[k],beta)*GThermal(r,t,beta,a);	
	}
}

// GaussianThermal2Disjoint
template <uint Dim>
void GaussianThermal2Disjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeighDisjoint(j,l.size());
		uint pk = posNeighDisjoint(k,l.size());
		number r = SpatialDistance(l[j],l[k]);
		number t = DXDisjoint(l,k,j,3,beta);
	
		result += f*(1.0+(number)(k<j))*DotDisjoint(l[pj],l[j],l[pk],l[k],beta)*HThermal(r,t,beta,a);	
	}
}

// GaussianThermalLRDisjoint
template <uint Dim>
void GaussianThermalLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		GaussianThermalDisjoint(j,k,l,beta,a,f,result);
}

// GaussianThermal2LRDisjoint
template <uint Dim>
void GaussianThermal2LRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		GaussianThermal2Disjoint(j,k,l,beta,a,f,result);
}

// InlineCurvatureMax
template <uint Dim>
void InlineCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	number temp = f*abs( Dot(l[pj],l[j],l[pj],l[j]) + Dot(l[nj],l[j],l[pj],l[j]) )/DistanceSquared(l[pj],l[j]);
	result = (temp>result? temp: result);
}

// InlineCurvatureMax
template <uint Dim>
void InlineCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = posNeigh(j,l.size());
		uint nj = negNeigh(j,l.size());
		Point<Dim> dd = l[pj] + l[nj] - 2.0*l[j];
		Point<Dim> d = l[pj] - l[j];
		number temp = f*abs(Dot(dd,d))/Dot(d,d);
		result = (temp>result? temp: result);
	}
}

// InlineCurvatureMaxDisjoint
template <uint Dim>
void InlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result){
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	number temp = f*abs( DotDisjoint(l[pj],l[j],l[pj],l[j],beta) + DotDisjoint(l[nj],l[j],l[pj],l[j],beta) )/DistanceSquaredDisjoint(l[pj],l[j],beta);
	result = (temp>result? temp: result);
}

// InlineCurvatureMaxDisjoint
template <uint Dim>
void InlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = posNeighDisjoint(j,l.size());
		uint nj = negNeighDisjoint(j,l.size());
		Point<Dim> dd = l[pj] + l[nj] - 2.0*l[j];
		Point<Dim> d = l[pj] - l[j];
		number temp = f*abs(DotDisjoint(dd,d,beta))/DotDisjoint(d,d,beta);
		result = (temp>result? temp: result);
	}
}

// CuspCurvatureMax
template <uint Dim>
void CuspCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	uint ppj = posNeigh(pj,l.size());
	uint nj = negNeigh(j,l.size());
	uint nnj = negNeigh(nj,l.size());
	Point<Dim> dd = l[ppj] + l[j] - 2.0*l[pj];
	Point<Dim> dp = l[pj] - l[j];
	Point<Dim> dn = l[j] - l[nj];
	number temp1 = f*abs(Dot(dd,dn)/Dot(dp,dn)), temp2;
	dd = l[nnj] + l[j] - 2.0*l[nj];
	temp2 = f*abs(Dot(dd,dp)/Dot(dp,dn));
	temp2 = (temp2>temp1? temp2: temp1);
	result = (temp2>result? temp2: result);
}

// CuspCurvatureMax
template <uint Dim>
void CuspCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = posNeigh(j,l.size());
		uint ppj = posNeigh(pj,l.size());
		uint nj = negNeigh(j,l.size());
		uint nnj = negNeigh(nj,l.size());
		Point<Dim> dd = l[ppj] + l[j] - 2.0*l[pj];
		Point<Dim> dp = l[pj] - l[j];
		Point<Dim> dn = l[j] - l[nj];
		number temp1 = f*abs(Dot(dd,dn)/Dot(dp,dn)), temp2;
		dd = l[nnj] + l[j] - 2.0*l[nj];
		temp2 = f*abs(Dot(dd,dp)/Dot(dp,dn));
		temp2 = (temp2>temp1? temp2: temp1);
		result = (temp2>result? temp2: result);
	}
}

// CuspCurvatureMaxDisjoint
template <uint Dim>
void CuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	uint ppj = posNeighDisjoint(pj,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	uint nnj = negNeighDisjoint(nj,l.size());
	Point<Dim> dd = l[ppj] + l[j] - 2.0*l[pj];
	Point<Dim> dp = l[pj] - l[j];
	Point<Dim> dn = l[j] - l[nj];
	number temp1 = f*abs(DotDisjoint(dd,dn,beta)/DotDisjoint(dp,dn,beta)), temp2;
	dd = l[nnj] + l[j] - 2.0*l[nj];
	temp2 = f*abs(DotDisjoint(dd,dp,beta)/DotDisjoint(dp,dn,beta));
	temp2 = (temp2>temp1? temp2: temp1);
	result = (temp2>result? temp2: result);
}

// CuspCurvatureMaxDisjoint
template <uint Dim>
void CuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = posNeighDisjoint(j,l.size());
		uint ppj = posNeighDisjoint(pj,l.size());
		uint nj = negNeighDisjoint(j,l.size());
		uint nnj = negNeighDisjoint(nj,l.size());
		Point<Dim> dd = l[ppj] + l[j] - 2.0*l[pj];
		Point<Dim> dp = l[pj] - l[j];
		Point<Dim> dn = l[j] - l[nj];
		number temp1 = f*abs(DotDisjoint(dd,dn,beta)/DotDisjoint(dp,dn,beta)), temp2;
		dd = l[nnj] + l[j] - 2.0*l[nj];
		temp2 = f*abs(DotDisjoint(dd,dp,beta)/DotDisjoint(dp,dn,beta));
		temp2 = (temp2>temp1? temp2: temp1);
		result = (temp2>result? temp2: result);
	}
}

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	Point<Dim> pp;
	number t_sqrd, t_dot_sqrd, t_t_dot, temp;
	uint pj = posNeigh(j,l.size());
	uint ppj = posNeigh(pj,l.size());
	pp = l[pj];
	pp *= 2.0;
	pp -= l[ppj];
	t_sqrd = DistanceSquared(l[pj],l[j]);
	t_dot_sqrd = DistanceSquared(l[j],pp);
	t_t_dot = Dot(l[j],pp,l[pj],l[j]);
	temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
	result =  (temp>result? temp: result);
}

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		Point<Dim> pp;
		number t_sqrd, t_dot_sqrd, t_t_dot, temp;
		uint pj = posNeigh(j,l.size());
		uint ppj = posNeigh(pj,l.size());
		pp = l[pj];
		pp *= 2.0;
		pp -= l[ppj];
		t_sqrd = DistanceSquared(l[pj],l[j]);
		t_dot_sqrd = DistanceSquared(l[j],pp);
		t_t_dot = Dot(l[j],pp,l[pj],l[j]);
		temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
		result = (temp>result? temp: result);
	}
}

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	number dx = DX(l,pj,j,Dim-2);
	number ddx = DX(l,pj,j,Dim-2) + DX(l,nj,j,Dim-2);
	number dy = DX(l,pj,j,Dim-1);
	number ddy = DX(l,pj,j,Dim-2) + DX(l,nj,j,Dim-1);
	number temp = f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5);
	result = (temp>result? temp: result);
}

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = posNeigh(j,l.size());
		uint nj = negNeigh(j,l.size());
		number dx = DX(l,pj,j,Dim-2);
		number ddx = DX(l,pj,j,Dim-2) + DX(l,nj,j,Dim-2);
		number dy = DX(l,pj,j,Dim-1);
		number ddy = DX(l,pj,j,Dim-2) + DX(l,nj,j,Dim-1);
		number temp = f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5);
		result = (temp>result? temp: result);
	}
}

// KGMaxDisjoint
template <uint Dim>
void KGMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	Point<Dim> pp;
	number t_sqrd, t_dot_sqrd, t_t_dot, temp;
	uint pj = posNeighDisjoint(j,l.size());
	uint ppj = posNeighDisjoint(pj,l.size());
	pp = l[pj];
	pp *= 2.0;
	pp -= l[ppj];
	t_sqrd = DistanceSquaredDisjoint(l[pj],l[j],beta);
	t_dot_sqrd = DistanceSquaredDisjoint(l[j],pp,beta);
	t_t_dot = DotDisjoint(l[j],pp,l[pj],l[j],beta);
	temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
	result =  (temp>result? temp: result);
}

// KGMaxDisjoint
template <uint Dim>
void KGMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		Point<Dim> pp;
		number t_sqrd, t_dot_sqrd, t_t_dot, temp;
		uint pj = posNeighDisjoint(j,l.size());
		uint ppj = posNeighDisjoint(pj,l.size());
		pp = l[pj];
		pp *= 2.0;
		pp -= l[ppj];
		t_sqrd = DistanceSquaredDisjoint(l[pj],l[j],beta);
		t_dot_sqrd = DistanceSquaredDisjoint(l[j],pp,beta);
		t_t_dot = DotDisjoint(l[j],pp,l[pj],l[j],beta);
		temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
		result = (temp>result? temp: result);
	}
}

// KGMaxPlaneDisjoint
template <uint Dim>
void KGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	number dx = DXDisjoint(l,pj,j,Dim-2,beta);
	number ddx = DXDisjoint(l,pj,j,Dim-2,beta) + DXDisjoint(l,nj,j,Dim-2,beta);
	number dy = DXDisjoint(l,pj,j,Dim-1,beta);
	number ddy = DXDisjoint(l,pj,j,Dim-2,beta) + DXDisjoint(l,nj,j,Dim-1,beta);
	number temp = f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5);
	result = (temp>result? temp: result);
}

// KGMaxPlaneDisjoint
template <uint Dim>
void KGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = posNeighDisjoint(j,l.size());
		uint nj = negNeighDisjoint(j,l.size());
		number dx = DXDisjoint(l,pj,j,Dim-2,beta);
		number ddx = DXDisjoint(l,pj,j,Dim-2,beta) + DXDisjoint(l,nj,j,Dim-2,beta);
		number dy = DXDisjoint(l,pj,j,Dim-1,beta);
		number ddy = DXDisjoint(l,pj,j,Dim-2,beta) + DXDisjoint(l,nj,j,Dim-1,beta);
		number temp = f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5);
		result = (temp>result? temp: result);
	}
}

// mdPX_nr
template<uint Dim>
void mdPX_nr(const Loop<Dim>& l,const uint& loc, const Point<Dim>& P, const number& f, vec& v) {
	for (uint mu=0; mu<Dim; mu++)
		v[loc*Dim+mu] += -f*P[mu];
}

// mdL_nr
template<uint Dim>
void mdL_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	
	number norm = Distance(l[j],l[nj]);
	number temp = DX(l,j,nj,mu)/norm;
	norm = Distance(l[pj],l[j]);
	temp += DX(l,j,pj,mu)/norm;
	
	v[j*Dim+mu] += -f*temp;
}

// mdLDisjoint_nr
template<uint Dim>
void mdLDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	
	number norm = DistanceDisjoint(l[j],l[nj],beta);
	number temp = DXDisjoint(l,j,nj,mu,beta)/norm;
	norm = DistanceDisjoint(l[pj],l[j],beta);
	temp += DXDisjoint(l,j,pj,mu,beta)/norm;
	
	v[j*Dim+mu] += -f*temp;
}

// ddL_nr
template<uint Dim>
void ddL_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	if (k==j) {
		number normn = Distance(l[j],l[nj]), normp = Distance(l[pj],l[j]), temp = 0.0;
		if (mu==nu)
			temp += 1.0/normn + 1.0/normp;
		temp -= DX(l,j,nj,mu)*DX(l,j,nj,nu)/pow(normn,3);
		temp -= DX(l,j,pj,mu)*DX(l,j,pj,nu)/pow(normp,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==nj) {
		number norm = Distance(l[j],l[nj]), temp = 0.0;
		if (mu==nu)
			temp -= 1.0/norm;
		temp += DX(l,j,nj,mu)*DX(l,j,nj,nu)/pow(norm,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==pj) {
		number norm = Distance(l[pj],l[j]), temp = 0.0;
		if (mu==nu)
			temp -= 1.0/norm;
		temp += DX(l,j,pj,mu)*DX(l,j,pj,nu)/pow(norm,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
}

// ddLDisjoint_nr
template<uint Dim>
void ddLDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
			const number& beta, const number& f, mat& m) {
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	if (k==j) {
		number normn = DistanceDisjoint(l[j],l[nj],beta), normp = DistanceDisjoint(l[pj],l[j],beta), temp = 0.0;
		if (mu==nu)
			temp += 1.0/normn + 1.0/normp;
		temp -= DXDisjoint(l,j,nj,mu,beta)*DXDisjoint(l,j,nj,nu,beta)/pow(normn,3);
		temp -= DXDisjoint(l,j,pj,mu,beta)*DXDisjoint(l,j,pj,nu,beta)/pow(normp,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==nj) {
		number norm = DistanceDisjoint(l[j],l[nj],beta), temp = 0.0;
		if (mu==nu)
			temp -= 1.0/norm;
		temp += DXDisjoint(l,j,nj,mu,beta)*DXDisjoint(l,j,nj,nu,beta)/pow(norm,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==pj) {
		number norm = DistanceDisjoint(l[pj],l[j],beta), temp = 0.0;
		if (mu==nu)
			temp -= 1.0/norm;
		temp += DXDisjoint(l,j,pj,mu,beta)*DXDisjoint(l,j,pj,nu,beta)/pow(norm,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
}

// mdDistPow_nr
template<uint Dim>
void mdDistPow_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& w,  const number& f, vec& v) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	
	number B_jmj = DistanceSquared(l[j],l[mj]);
	number B_jpj = DistanceSquared(l[j],l[pj]);
	
	number temp =  w*pow(B_jmj,-1 + w/2.)*DX(l,j,mj,mu) \
				+ w*pow(B_jpj,-1 + w/2.)*DX(l,j,pj,mu);
	
	v[j*Dim+mu] += -f*temp;
}

// ddDistPow_nr
template<uint Dim>
void ddDistPow_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& w\
				, const number& f, mat& m) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	
	number B_jmj, B_jpj, temp = 0.0; 
	
	if (k==j) {
		B_jmj = DistanceSquared(l[j],l[mj]);
		B_jpj = DistanceSquared(l[j],l[pj]);
		if (mu==nu)
			temp += w*pow(B_jmj,-1 + w/2.) \
 					+ w*pow(B_jpj,-1 + w/2.);
		temp += - 2.0*w*pow(B_jmj,-2 + w/2.)*DX(l,j,mj,mu)*DX(l,j,mj,nu) \
				 + pow(w,2)*pow(B_jmj,-2 + w/2.)*DX(l,j,mj,mu)*DX(l,j,mj,nu) \
				 - 2.0*w*pow(B_jpj,-2 + w/2.)*DX(l,j,pj,mu)*DX(l,j,pj,nu) \
				 + pow(w,2)*pow(B_jpj,-2 + w/2.)*DX(l,j,pj,mu)*DX(l,j,pj,nu);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==mj) {
		B_jmj = DistanceSquared(l[j],l[mj]);
		if (mu==nu)
			temp += - w*pow(B_jmj,-1 + w/2.);
		temp += 2.0*w*pow(B_jmj,-2 + w/2.)*DX(l,j,mj,mu)*DX(l,j,mj,nu) \
 				- pow(w,2)*pow(B_jmj,-2 + w/2.)*DX(l,j,mj,mu)*DX(l,j,mj,nu);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==pj) {
		B_jpj = DistanceSquared(l[j],l[pj]);
		if (mu==nu)
			temp +=  - w*pow(B_jpj,-1 + w/2.);
		temp += + 2.0*w*pow(B_jpj,-2 + w/2.)*DX(l,j,pj,mu)*DX(l,j,pj,nu) \
 				- pow(w,2)*pow(B_jpj,-2 + w/2.)*DX(l,j,pj,mu)*DX(l,j,pj,nu);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
}

// mdI0_nr
template<uint Dim>
void mdI0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
	if (mu==(Dim-1)) {
		uint nj = negNeigh(j,l.size());
		v[j*Dim+mu] += -f*DX(l,nj,j,Dim-2);
	}
	else if (mu==(Dim-2)) {
		uint pj = posNeigh(j,l.size());
		v[j*Dim+mu] += -f*DX(l,pj,j,Dim-1);
	}
}

// mdI0Disjoint_nr
template<uint Dim>
void mdI0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	if (mu==(Dim-1)) {
		uint nj = negNeighDisjoint(j,l.size());
		v[j*Dim+mu] += -f*DXDisjoint(l,nj,j,Dim-2,beta);
	}
	else if (mu==(Dim-2)) {
		uint pj = posNeighDisjoint(j,l.size());
		v[j*Dim+mu] += -f*DXDisjoint(l,pj,j,Dim-1,beta);
	}
}

// ddI0_nr
template<uint Dim>
void ddI0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
	if (mu==(Dim-1) && nu==(Dim-2)) {
		uint nj = negNeigh(j,l.size());
		if (k==nj)
			m(j*Dim+mu,k*Dim+nu) += f;
		else if (k==j)
			m(j*Dim+mu,k*Dim+nu) -= f;
	}
	else if (mu==(Dim-2) && nu==(Dim-1)) {
		uint pj = posNeigh(j,l.size());
		if (k==pj)
			m(j*Dim+mu,k*Dim+nu) += f;
		else if (k==j)
			m(j*Dim+mu,k*Dim+nu) -= f;
	}
}

// ddI0Disjoint_nr
template<uint Dim>
void ddI0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& beta, const number& f, mat& m) {
	if (mu==(Dim-1) && nu==(Dim-2)) {
		uint nj = negNeighDisjoint(j,l.size());
		if (k==nj)
			m(j*Dim+mu,k*Dim+nu) += f;
		else if (k==j)
			m(j*Dim+mu,k*Dim+nu) -= f;
	}
	else if (mu==(Dim-2) && nu==(Dim-1)) {
		uint pj = posNeighDisjoint(j,l.size());
		if (k==pj)
			m(j*Dim+mu,k*Dim+nu) += f;
		else if (k==j)
			m(j*Dim+mu,k*Dim+nu) -= f;
	}
}


// mdIn_nr
template<uint Dim>
void mdIn_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& f, vec& v) {
	if (abs(n)>MIN_NUMBER) {
		if (mu==(Dim-1)) {
			uint nj = negNeigh(j,l.size());
			v[j*Dim+mu] += -f*(pow((l[nj])[Dim-2],n)-pow((l[j])[Dim-2],n));
		}
		else if (mu==(Dim-2)) {
			uint pj = posNeigh(j,l.size());
			v[j*Dim+mu] += -f*n*pow((l[j])[Dim-2],n-1.0)*DX(l,pj,j,Dim-1);
		}
	}
}

// mdInDisjoint_nr
template<uint Dim>
void mdInDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, vec& v) {
	if (abs(n)>MIN_NUMBER) {
		if (mu==(Dim-1)) {
			uint nj = negNeighDisjoint(j,l.size());
			v[j*Dim+mu] += -f*(pow((l[nj])[Dim-2],n)-pow((l[j])[Dim-2],n));
		}
		else if (mu==(Dim-2)) {
			uint pj = posNeighDisjoint(j,l.size());
			v[j*Dim+mu] += -f*n*pow((l[j])[Dim-2],n-1.0)*DXDisjoint(l,pj,j,Dim-1,beta);
		}
	}
}

// ddIn_nr
template<uint Dim>
void ddIn_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& f, mat& m) {
	if (abs(n)>MIN_NUMBER) {
		if (mu==(Dim-2) && nu==(Dim-2) && abs(n-1.0)>MIN_NUMBER) {
			uint pj = posNeigh(j,l.size());
			if (k==j)
				m(j*Dim+mu,k*Dim+nu) += f*n*(n-1.0)*pow((l[j])[Dim-2],n-2.0)*DX(l,pj,j,Dim-1);
		}
		else if (mu==(Dim-1) && nu==(Dim-2)) {
			uint nj = negNeigh(j,l.size());
			if (k==nj)
				m(j*Dim+mu,k*Dim+nu) += f*n*pow((l[k])[Dim-2],n-1.0);
			else if (k==j)
				m(j*Dim+mu,k*Dim+nu) -= f*n*pow((l[k])[Dim-2],n-1.0);
		}
		else if (mu==(Dim-2) && nu==(Dim-1)) {
			uint pj = posNeigh(j,l.size());
			if (k==pj)
				m(j*Dim+mu,k*Dim+nu) += f*n*pow((l[j])[Dim-2],n-1.0);
			else if (k==j)
				m(j*Dim+mu,k*Dim+nu) -= f*n*pow((l[j])[Dim-2],n-1.0);
		}
	}
}

// ddInDisjoint_nr
template<uint Dim>
void ddInDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, mat& m) {
	if (abs(n)>MIN_NUMBER) {
		if (mu==(Dim-2) && nu==(Dim-2) && abs(n-1.0)>MIN_NUMBER) {
			uint pj = posNeighDisjoint(j,l.size());
			if (k==j)
				m(j*Dim+mu,k*Dim+nu) += f*n*(n-1.0)*pow((l[j])[Dim-2],n-2.0)*DXDisjoint(l,pj,j,Dim-1,beta);
		}
		else if (mu==(Dim-1) && nu==(Dim-2)) {
			uint nj = negNeighDisjoint(j,l.size());
			if (k==nj)
				m(j*Dim+mu,k*Dim+nu) += f*n*pow((l[k])[Dim-2],n-1.0);
			else if (k==j)
				m(j*Dim+mu,k*Dim+nu) -= f*n*pow((l[k])[Dim-2],n-1.0);
		}
		else if (mu==(Dim-2) && nu==(Dim-1)) {
			uint pj = posNeighDisjoint(j,l.size());
			if (k==pj)
				m(j*Dim+mu,k*Dim+nu) += f*n*pow((l[j])[Dim-2],n-1.0);
			else if (k==j)
				m(j*Dim+mu,k*Dim+nu) -= f*n*pow((l[j])[Dim-2],n-1.0);
		}
	}
}

// mdS0_nr
template<uint Dim>
void mdS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	v[j*Dim+mu] += -f*(DX(l,j,nj,mu) + DX(l,j,pj,mu))*(number)l.size()/2.0;
}

// mdS0Disjoint_nr
template<uint Dim>
void mdS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	v[j*Dim+mu] += -f*(DXDisjoint(l,j,nj,mu,beta) + DXDisjoint(l,j,pj,mu,beta))*(number)l.size()/2.0;
}

// ddS0_nr
template<uint Dim>
void ddS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
	if (mu==nu) {
		uint pj = posNeigh(j,l.size());
		uint nj = negNeigh(j,l.size());
		if (k==j) {
			m(j*Dim+mu,k*Dim+nu) += f*(number)l.size();
		}
		else if (k==nj) {
			m(j*Dim+mu,k*Dim+nu) -= f*0.5*(number)l.size();
		}
		else if (k==pj) {
			m(j*Dim+mu,k*Dim+nu) -= f*0.5*(number)l.size();
		}
	}
}

// ddS0Disjoint_nr
template<uint Dim>
void ddS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
	if (mu==nu) {
		uint pj = posNeighDisjoint(j,l.size());
		uint nj = negNeighDisjoint(j,l.size());
		if (k==j) {
			m(j*Dim+mu,k*Dim+nu) += f*(number)l.size();
		}
		else if (k==nj) {
			m(j*Dim+mu,k*Dim+nu) -= f*0.5*(number)l.size();
		}
		else if (k==pj) {
			m(j*Dim+mu,k*Dim+nu) -= f*0.5*(number)l.size();
		}
	}
}

// mdsqrtS0_nr
template<uint Dim>
void mdsqrtS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0, const number& f, vec& v) {
	mdS0_nr(j,mu,l,f*2.0/sqrt4s0,v);
}

// mdsqrtS0Disjoint_nr
template<uint Dim>
void mdsqrtS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0,\
 const number& beta, const number& f, vec& v) {
 	mdS0Disjoint_nr(j,mu,l,beta,f*2.0/sqrt4s0,v);
 }

// ddsqrtS0_nr
template<uint Dim>
void ddsqrtS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& f, mat& m) {
	ddS0_nr(j,mu,k,nu,l,f*2.0/sqrt4s0,m);
	
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	uint pk = posNeigh(k,l.size());
	uint nk = negNeigh(k,l.size());
	
	m(j*Dim+mu,k*Dim+nu) -= (f*pow((number)l.size(),2)/pow(sqrt4s0,3)) * (DX(l,j,nj,mu) + DX(l,j,pj,mu)) \
												* (DX(l,k,nk,nu) + DX(l,k,pk,nu));
}


// ddsqrtS0Disjoint_nr
template<uint Dim>
void ddsqrtS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& beta, const number& f, mat& m) {
	ddS0Disjoint_nr(j,mu,k,nu,l,f*2.0/sqrt4s0,m);
	
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	uint nk = negNeighDisjoint(k,l.size());
	
	m(j*Dim+mu,k*Dim+nu) -= (f*pow((number)l.size(),2)/pow(sqrt4s0,3)) * (DXDisjoint(l,j,nj,mu,beta) + DXDisjoint(l,j,pj,mu,beta)) \
												* (DXDisjoint(l,k,nk,nu,beta) + DXDisjoint(l,k,pk,nu,beta));				
}

// mdVor_nr
template<uint Dim>
void mdVor_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdVor_nr( j, mu, i, l, a, f, v);
	
}

// mdVor_nr
template<uint Dim>
void mdVor_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	cerr << "mdVor_nr Error: no script written for dim = " << Dim << endl;
}

// ddVor_nr
template<uint Dim>
void ddVor_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "mddVor_nr Error: no script written for dim = " << Dim << endl;
}

// mdVlr_nr
template<uint Dim>
void mdVlr_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdVlr_nr( j, mu, i, l, a, f, v);
	
}

// mdVlr_nr
template<uint Dim>
void mdVlr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	cerr << "mdVlr_nr Error: no script written for dim = " << Dim << endl;
}

// ddVlr_nr
template<uint Dim>
void ddVlr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "mddVlr_nr Error: no script written for dim = " << Dim << endl;
}

// mdVer_nr
template<uint Dim>
void mdVer_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdVer_nr( j, mu, i, l, a, f, v);
	
}

// mdVer_nr
template<uint Dim>
void mdVer_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	cerr << "mdVer_nr Error: no script written for dim = " << Dim << endl;
}

// ddVer_nr
template<uint Dim>
void ddVer_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "mddVer_nr Error: no script written for dim = " << Dim << endl;
}

// mdVdr_nr
template<uint Dim>
void mdVdr_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdVdr_nr( j, mu, i, l, a, f, v);
	
}

// mdVdr_nr
template<uint Dim>
void mdVdr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	cerr << "mdVdr_nr Error: no script written for dim = " << Dim << endl;
}

// ddVdr_nr
template<uint Dim>
void ddVdr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "mddVdr_nr Error: no script written for dim = " << Dim << endl;
}

// mdVthr_nr
template<uint Dim>
void mdVthr_nr(const uint& j, const uint& mu, const Loop<Dim>& l,\
		 const number& theta, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdVthr_nr( j, mu, i, l, a, f, v);
	
}

// mdVthr_nr
template<uint Dim>
void mdVdr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
		 const number& theta, const number& a, const number& f, vec& v) {
	cerr << "mdVthr_nr Error: no script written for dim = " << Dim << endl;
}

// mdVnonrelDisjoint_nr
template<uint Dim>
void mdVnonrelDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	cerr << "mdVnonrelDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// mdVnonrelrDisjoint_nr
template<uint Dim>
void mdVnonrelrDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& a\
						, const number& f, vec& v) {
	cerr << "mdVnonrelrDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// mdVthr_nr
template<uint Dim>
void mdVthr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mddVthr_nr Error: no script written for dim = " << Dim << endl;
}

// ddVthr_nr
template<uint Dim>
void ddVthr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& theta, const number& a, const number& f, mat& m) {
	cerr << "mddVthr_nr Error: no script written for dim = " << Dim << endl;
}

// mdVthrDisjoint_nr
template<uint Dim>
void mdVthrDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mddVthrDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// mdVthrDisjointLR_nr
template<uint Dim>
void mdVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mddVthrDisjointLR_nr Error: no script written for dim = " << Dim << endl;
}

// ddVthrDisjoint_nr
template<uint Dim>
void ddVthrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddVthrDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// ddVthrDisjointLR_nr
template<uint Dim>
void ddVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddVthrDisjointLR_nr Error: no script written for dim = " << Dim << endl;
}

// ddVnonrelDisjoint_nr
template<uint Dim>
void ddVnonrelDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& f, mat& m) {
	cerr << "ddVnonrelDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussian_nr
template<uint Dim>
void mdGaussian_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdGaussian_nr( j, mu, i, l, a, f, v);
	
}

// mdGaussian_nr
template<uint Dim>
void mdGaussian_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	cerr << "mdGaussian_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussianThermal_nr
template<uint Dim>
void mdGaussianThermal_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianThermal_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussianThermal2_nr
template<uint Dim>
void mdGaussianThermal2_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianThermal2_nr Error: no script written for dim = " << Dim << endl;
}

// ddGaussian_nr
template<uint Dim>
void ddGaussian_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "ddGaussian_nr Error: no script written for dim = " << Dim << endl;
}


// ddGaussianThermal_nr
template<uint Dim>
void ddGaussianThermal_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddGaussianThermal_nr Error: no script written for dim = " << Dim << endl;
}

// ddGaussianThermal2_nr
template<uint Dim>
void ddGaussianThermal2_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddGaussianThermal2_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussianDisjoint_nr
template<uint Dim>
void mdGaussianDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	for (uint i=0; i<l.size(); i++)
		mdGaussianDisjoint_nr( j, mu, i, l, beta, a, f, v);
}

// mdGaussianDisjoint_nr
template<uint Dim>
void mdGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussianLRDisjoint_nr
template<uint Dim>
void mdGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianLRDisjoint_nr Error: no script written for dim = " << Dim << endl;
}


// mdGaussianThermalDisjoint_nr
template<uint Dim>
void mdGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianThermalDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussianThermal2Disjoint_nr
template<uint Dim>
void mdGaussianThermal2Disjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianThermal2Disjoint_nr Error: no script written for dim = " << Dim << endl;
}
			
// mdGaussianThermalLRDisjoint_nr
template<uint Dim>
void mdGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianThermalLRDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussianThermal2LRDisjoint_nr
template<uint Dim>
void mdGaussianThermal2LRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianThermal2LRDisjoint_nr Error: no script written for dim = " << Dim << endl;
}
		 
// ddGaussianDisjoint_nr
template<uint Dim>
void ddGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddGaussianDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// ddGaussianLRDisjoint_nr
template<uint Dim>
void ddGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddGaussianLRDisjoint_nr Error: no script written for dim = " << Dim << endl;
}
						 
// ddGaussianThermalDisjoint_nr
template<uint Dim>
void ddGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddGaussianThermalDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

						 
// ddGaussianThermalLRDisjoint_nr
template<uint Dim>
void ddGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddGaussianThermalLRDisjoint_nr Error: no script written for dim = " << Dim << endl;	 
}

// mdGamma_nr
template<uint Dim>
void mdFGamma_nr(const Loop<Dim>& l, const uint& j, const number& f, vec& v) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	
	Point<Dim> P = l[j]-l[mj];
	Point<Dim> Q = l[pj]-l[j];
	
	number p = Norm(P);
	number q = Norm(Q);
	number pq = Dot(P,Q);
	number c = pq/p/q;
	number s = sqrt(1.0-c*c);
	number DG = -c/pow(s,2) + acos(c)/pow(s,3);
	//number dcdp, dcdq;
	
	for (uint mu=0; mu<Dim; mu++) {
		/*dcdp = Q[mu]/p/q - pq*P[mu]/q/pow(p,3);
		dcdq = P[mu]/p/q - pq*Q[mu]/p/pow(q,3);
		
		v(mj*Dim+mu) += -f*(-DG*dcdp);
		v(j*Dim+mu) += -f*(DG*(dcdp - dcdq));
		v(pj*Dim+mu) += -f*(DG*dcdq);*/
		v(mj*Dim+mu) += -f*( (DG*(pq*P[mu] - pow(p,2)*Q[mu]))/(pow(p,3)*q) );
		v(j*Dim+mu) += -f*( -((DG*((pow(p,2) + pq)*pow(q,2)*P[mu] - pow(p,2)*(pq + pow(q,2))*Q[mu]))/(pow(p,3)*pow(q,3))) );
		v(pj*Dim+mu) += -f*( (DG*(pow(q,2)*P[mu] - pq*Q[mu]))/(p*pow(q,3)) );
	}
}

// ddGamma_nr
template<uint Dim>
void ddFGamma_nr(const Loop<Dim>& l, const uint& j, const number& f, mat& m) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	
	Point<Dim> P = l[j]-l[mj];
	Point<Dim> Q = l[pj]-l[j];
	
	number p = Norm(P);
	number q = Norm(Q);
	number pq = Dot(P,Q);
	number c = pq/p/q;
	number s = sqrt(1.0-c*c);
	number DG = -c/pow(s,2) + acos(c)/pow(s,3);
	number DDG = (-s*(2.0+c*c) + 3.0*c*acos(c))/pow(s,5);
	//number dcdp_mu, dcdp_nu, dcdq_mu, dcdq_nu, ddcdpdp, ddcdpdq, ddcdqdp, ddcdqdq; // mu before nu in ordering
	
	for (uint mu=0; mu<Dim; mu++) {
		/*dcdp_mu = Q[mu]/p/q - pq*P[mu]/q/pow(p,3);
		dcdq_mu = P[mu]/p/q - pq*Q[mu]/p/pow(q,3);*/
		
		for (uint nu=0; nu<Dim; nu++) {
			/*dcdp_nu = Q[nu]/p/q - pq*P[nu]/q/pow(p,3);
			dcdq_nu = P[nu]/p/q - pq*Q[nu]/p/pow(q,3);
			
			ddcdpdp = -(Q[mu]*P[nu]+P[mu]*Q[nu])/pow(p,3)/q + 3.0*pq*P[mu]*P[nu]/pow(p,5)/q;
			ddcdpdp += (mu==nu? -pq/pow(p,3)/q: 0.0);
			
			ddcdqdq = -(Q[mu]*P[nu]+P[mu]*Q[nu])/p/pow(q,3) + 3.0*pq*Q[mu]*Q[nu]/p/pow(q,5);
			ddcdqdq += (mu==nu? -pq/p/pow(q,3): 0.0);
			
			ddcdpdq = -Q[mu]*Q[nu]/p/pow(q,3) - P[mu]*P[mu]/pow(p,3)/q + pq*P[mu]*Q[nu]/pow(p,3)/pow(q,3);
			ddcdpdq += (mu==nu? 1.0/p/q: 0.0);
			
			ddcdqdp = -Q[mu]*Q[nu]/p/pow(q,3) - P[mu]*P[nu]/pow(p,3)/q + pq*Q[mu]*P[nu]/pow(p,3)/pow(q,3);
			ddcdqdp += (mu==nu? 1.0/p/q: 0.0);
			
			
			m(mj*Dim+mu,mj*Dim+nu) += f*( DDG*dcdp_mu*dcdp_nu \
										+ DG*ddcdpdp );	
			m(mj*Dim+mu,j*Dim+nu) += f*( -DDG*dcdp_mu*(dcdp_nu-dcdq_nu) \
										- DG*(ddcdpdp - ddcdpdq) );						
			m(mj*Dim+mu,pj*Dim+nu) += f*( -DDG*dcdp_mu*dcdq_nu \
										- DG*ddcdpdq );		
			m(j*Dim+mu,mj*Dim+nu) += f*( -DDG*(dcdp_mu-dcdq_mu)*dcdp_nu \
									- DG*(ddcdpdp - ddcdqdp) );
			m(j*Dim+mu,j*Dim+nu) += f*( DDG*(dcdp_mu-dcdq_mu)*(dcdp_nu-dcdq_nu) \
									+ DG*(ddcdpdp - ddcdpdq - ddcdqdp + ddcdqdq) );
			m(j*Dim+mu,pj*Dim+nu) += f*( DDG*(dcdp_mu-dcdq_mu)*dcdq_nu \
									- DG*(ddcdpdq + ddcdqdq) );
			m(pj*Dim+mu,mj*Dim+nu) += f*( -DDG*dcdq_mu*dcdp_nu \
										- DG*ddcdqdp );		
			m(pj*Dim+mu,j*Dim+nu) += f*( DDG*dcdq_mu*(dcdp_nu-dcdq_nu) \
										- DG*(-ddcdqdp + ddcdqdq) );		
			m(pj*Dim+mu,pj*Dim+nu) += f*( DDG*dcdq_mu*dcdq_nu \
										+ DG*ddcdqdq );	*/
			m(mj*Dim+mu,mj*Dim+nu) += f*( (P[mu]*(pq*(DDG*pq + 3.0*DG*p*q)*P[nu] - pow(p,2)*(DDG*pq + DG*p*q)*Q[nu]) + \
     pow(p,2)*(-((DDG*pq + DG*p*q)*P[nu]*Q[mu]) + p*(DDG*p*Q[mu]*Q[nu] - DG*pq*q*delta(mu,nu))))/(pow(p,6)*pow(q,2)) );	
     
			m(mj*Dim+mu,j*Dim+nu) += f*( (P[mu]*(-(pow(q,2)*(DDG*pq*(pow(p,2) + pq) + DG*p*(pow(p,2) + 3.0*pq)*q)*P[nu]) + \
        pow(p,2)*(DDG*pq + DG*p*q)*(pq + pow(q,2))*Q[nu]) + \
     pow(p,2)*(pow(q,2)*(DDG*(pow(p,2) + pq) + DG*p*q)*P[nu]*Q[mu] - \
        p*(p*(DG*p*q + DDG*(pq + pow(q,2)))*Q[mu]*Q[nu] - DG*(pow(p,2) + pq)*pow(q,3)*delta(mu,nu))))/ \
   (pow(p,6)*pow(q,4)) );			
   			
			m(mj*Dim+mu,pj*Dim+nu) += f*( ((-(DDG*pq) - DG*p*q)*P[mu]*(-(pow(q,2)*P[nu]) + pq*Q[nu]) + \
     pow(p,2)*(-(DDG*pow(q,2)*P[nu]*Q[mu]) + (DDG*pq + DG*p*q)*Q[mu]*Q[nu] - DG*p*pow(q,3)*delta(mu,nu)))/ \
   (pow(p,4)*pow(q,4)) );	
											
			m(j*Dim+mu,mj*Dim+nu) += f*( (pow(q,2)*P[mu]*(-((DDG*pq*(pow(p,2) + pq) + DG*p*(pow(p,2) + 3.0*pq)*q)*P[nu]) + \
        pow(p,2)*(DDG*(pow(p,2) + pq) + DG*p*q)*Q[nu]) + \
     pow(p,2)*((DDG*pq + DG*p*q)*(pq + pow(q,2))*P[nu]*Q[mu] - \
        p*(p*(DG*p*q + DDG*(pq + pow(q,2)))*Q[mu]*Q[nu] - DG*(pow(p,2) + pq)*pow(q,3)*delta(mu,nu))))/ \
   (pow(p,6)*pow(q,4)) );
									
			m(j*Dim+mu,j*Dim+nu) += f*( (pow(q,2)*P[mu]*(pow(q,2)*(DDG*pow(pow(p,2) + pq,2) \
			+ DG*p*(2.0*pow(p,2) + 3.0*pq)*q)*P[nu] - 
        pow(p,2)*(DDG*(pow(p,2) + pq)*(pq + pow(q,2)) + DG*p*q*(pow(p,2) + pq + pow(q,2)))*Q[nu]) - \
     pow(p,2)*(pow(q,2)*(DDG*(pow(p,2) + pq)*(pq + pow(q,2)) + DG*p*q*(pow(p,2) + pq + pow(q,2)))*P[nu]* \
         Q[mu] - p*(p*(DDG*pow(pq + pow(q,2),2) + DG*p*q*(3.0*pq + 2.0*pow(q,2)))*Q[mu]*Q[nu] - \
           DG*pow(q,3)*(pq*pow(q,2) + pow(p,2)*(pq + 2.0*pow(q,2)))*delta(mu,nu))))/(pow(p,6)*pow(q,6)) );
									
			m(j*Dim+mu,pj*Dim+nu) += f*( (pow(q,2)*P[mu]*(-(pow(q,2)*(DDG*(pow(p,2) + pq) + DG*p*q)*P[nu]) \
			+ (pow(p,2) + pq)*(DDG*pq + DG*p*q)*Q[nu]) + \
     pow(p,2)*(pow(q,2)*(DG*p*q + DDG*(pq + pow(q,2)))*P[nu]*Q[mu] - \
        (DDG*pq*(pq + pow(q,2)) + DG*p*q*(3.0*pq + pow(q,2)))*Q[mu]*Q[nu]\
         + DG*p*pow(q,3)*(pq + pow(q,2))*delta(mu,nu)))/ \
   (pow(p,4)*pow(q,6)) );
									
			m(pj*Dim+mu,mj*Dim+nu) += f*( (-(pq*(DDG*pq + DG*p*q)*P[nu]*Q[mu]) \
			+ pow(q,2)*P[mu]*((DDG*pq + DG*p*q)*P[nu] - DDG*pow(p,2)*Q[nu]) + \
     pow(p,2)*((DDG*pq + DG*p*q)*Q[mu]*Q[nu] - DG*p*pow(q,3)*delta(mu,nu)))/(pow(p,4)*pow(q,4)) );		
										
			m(pj*Dim+mu,j*Dim+nu) += f*( ((pow(p,2) + pq)*pow(q,2)*(DDG*pq + DG*p*q)*P[nu]*Q[mu] + \
     pow(q,2)*P[mu]*(-(pow(q,2)*(DDG*(pow(p,2) + pq) + DG*p*q)*P[nu]) + \
        pow(p,2)*(DG*p*q + DDG*(pq + pow(q,2)))*Q[nu]) - \
     pow(p,2)*((DDG*pq*(pq + pow(q,2)) + DG*p*q*(3.0*pq + pow(q,2)))*Q[mu]*Q[nu] - \
        DG*p*pow(q,3)*(pq + pow(q,2))*delta(mu,nu)))/(pow(p,4)*pow(q,6)) );	
											
			m(pj*Dim+mu,pj*Dim+nu) += f*( (-(pow(q,2)*(DDG*pq + DG*p*q)*P[nu]*Q[mu]) \
			+ pow(q,2)*P[mu]*(DDG*pow(q,2)*P[nu] - (DDG*pq + DG*p*q)*Q[nu]) + \
     pq*((DDG*pq + 3.0*DG*p*q)*Q[mu]*Q[nu] - DG*p*pow(q,3)*delta(mu,nu)))/(pow(p,2)*pow(q,6)) );							
		}
	}
}

// ErgS0_nr
template<uint Dim>
void ErgS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, number& erg) {
	erg += f*DX(l,j,mu)*(number)l.size()/2.0;
}

// ErgS0Disjoint_nr
template<uint Dim>
void ErgS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, number& erg) {
	erg += f*DXDisjoint(l,j,mu,beta)*(number)l.size()/2.0;
}

// ErgL_nr
template<uint Dim>
void ErgL_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, number& erg) {
	uint pj = posNeigh(j,l.size());
	erg += f*DX(l,pj,j,mu)/Distance(l[pj],l[j]);
}

// ErgLDisjoint_nr
template<uint Dim>
void ErgLDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, number& erg) {
	uint pj = posNeighDisjoint(j,l.size());
	erg += f*DXDisjoint(l,pj,j,mu,beta)/DistanceDisjoint(l[pj],l[j],beta);
}

// ErgsqrtS0_nr
template<uint Dim>
void ErgsqrtS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& f, number& erg) {
	ErgS0_nr(l,j,mu,f*2.0/sqrt4s0,erg);
}

// ErgsqrtS0Disjoint_nr
template<uint Dim>
void ErgsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, number& erg) {
	ErgS0Disjoint_nr(l,j,mu,beta,f*2.0/sqrt4s0,erg);
}

// ErgI0_nr
template<uint Dim>
void ErgI0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeigh(j,l.size());
		erg += -f*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*(l[j])[Dim-2];
}

// ErgI0Disjoint_nr
template<uint Dim>
void ErgI0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeighDisjoint(j,l.size());
		erg += -f*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*(l[j])[Dim-2];
}

// ErgIn_nr
template<uint Dim>
void ErgIn_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeigh(j,l.size());
		erg += -f*n*pow((l[j])[Dim-2],n-1.0)*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*pow((l[j])[Dim-2],n);
}

// ErgInDisjoint_nr
template<uint Dim>
void ErgInDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& beta, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeighDisjoint(j,l.size());
		erg += -f*n*pow((l[j])[Dim-2],n-1.0)*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*pow((l[j])[Dim-2],n);
}

// ErgVor_nr
template<uint Dim>
void ErgVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg) {
	erg += f*DX(l,k,mu)/(DistanceSquared(l[j],l[k]) + a*a);
}

// ErgVthr_nr
template<uint Dim>
void ErgVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);

	erg += f*(-pow(2.0*PI,2))*DX(l,k,mu)*FThermal(r,t,beta,a);
}

// ErgVthrDisjoint_nr
template<uint Dim>
void ErgVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1); // intentionally not DXDisjoint(l,k,j,Dim-1,beta);

	erg += f*(-pow(2.0*PI,2))*DXDisjoint(l,k,mu,beta)*FThermal(r,t,beta,a);	
}

// ErgVthrDisjointLR_nr
template<uint Dim>
void ErgVthrDisjointLR_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		ErgVthrDisjoint_nr(l,j,mu,k,beta,a,f,erg);
	else
		ErgVthrDisjoint_nr(l,j,mu,k,beta,0.0,f,erg);
}

// ErgGaussian_nr
template<uint Dim>
void ErgGaussian_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg) {
	erg += f*DX(l,k,mu)*exp(-DistanceSquared(l[j],l[k])/a/a);
}

// ErgGaussianDisjoint_nr
template<uint Dim>
void ErgGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	erg += f*DXDisjoint(l,k,mu,beta)*exp(-DistanceSquared(l[j],l[k])/a/a); //intentionally not disjoint
}

// ErgGaussianLRDisjoint_nr
template<uint Dim>
void ErgGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		ErgGaussianDisjoint_nr(l,j,mu,k,beta,a,f,erg);
}

// ErgGaussianThermal_nr
template<uint Dim>
void ErgGaussianThermal_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DX(l,k,mu)*GThermal(r,t,beta,a);
}

// ErgGaussianThermal2_nr
template<uint Dim>
void ErgGaussianThermal2_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DX(l,k,mu)*HThermal(r,t,beta,a);
}

// ErgGaussianThermalDisjoint_nr
template<uint Dim>
void ErgGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DXDisjoint(l,k,mu,beta)*GThermal(r,t,beta,a);
}

// ErgGaussianThermal2Disjoint_nr
template<uint Dim>
void ErgGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DXDisjoint(l,k,mu,beta)*HThermal(r,t,beta,a);
}

// ErgGaussianThermalLRDisjoint_nr
template<uint Dim>
void ErgGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		ErgGaussianThermalDisjoint_nr(l,j,mu,k,beta,a,f,erg);
}

// ErgGaussianThermal2LRDisjoint_nr
template<uint Dim>
void ErgGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		ErgGaussianThermal2Disjoint_nr(l,j,mu,k,beta,a,f,erg);
}

// PS0_nr
template<uint Dim>
void PS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, vec& v) {
	ErgS0_nr(l,j,mu,f,v[Dim*j+mu]);
}

// PS0Disjoint_nr
template<uint Dim>
void PS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, vec& v) {
	ErgS0Disjoint_nr(l,j,mu,beta,f,v[Dim*j+mu]);
}

// PL_nr
template<uint Dim>
void PL_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, vec& v) {
	ErgL_nr(l,j,mu,f,v[Dim*j+mu]);
}

// PLDisjoint_nr
template<uint Dim>
void PLDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, vec& v) {
	ErgLDisjoint_nr(l,j,mu,beta,f,v[Dim*j+mu]);
}

// PsqrtS0_nr
template<uint Dim>
void PsqrtS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& f, vec& v) {
	PS0_nr(l,j,mu,f*2.0/sqrt4s0,v);
}

// PsqrtS0Disjoint_nr
template<uint Dim>
void PsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, vec& v) {
	PS0Disjoint_nr(l,j,mu,beta,f*2.0/sqrt4s0,v);
}

// PI0_nr
template<uint Dim>
void PI0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, vec& v) {
	ErgI0_nr(l,j,mu,f,v[Dim*j+mu]);
}

// PI0Disjoint_nr
template<uint Dim>
void PI0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, vec& v) {
	ErgI0Disjoint_nr(l,j,mu,beta,f,v[Dim*j+mu]);
}

// PIn_nr
template<uint Dim>
void PIn_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& f, vec& v) {
	ErgIn_nr(l,j,mu,n,f,v[Dim*j+mu]);
}

// PInDisjoint_nr
template<uint Dim>
void PInDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& beta, const number& f, vec& v) {
	ErgInDisjoint_nr(l,j,mu,n,beta,f,v[Dim*j+mu]);
}

// PVor_nr
template<uint Dim>
void PVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v) {
	ErgVor_nr(l,j,mu,k,a,2.0*f,v[Dim*j+mu]);
}

// PRVor_nr
template<uint Dim>
void PRVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v) {
	number res = 0.0;
	for (uint i=0; i<=j; i++) {
		res += DX(l,k,nu)*DX(l,i,nu)*(-2.0*DX(l,k,i,mu))/pow(DistanceSquared(l[i],l[k]) + a*a,2);		
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PVthr_nr
template<uint Dim>
void PVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {	
	ErgVthr_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRVthr_nr
template<uint Dim>
void PRVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeigh(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeigh(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += (-pow(2.0*PI,2))*DX(l,pk,k,nu)*DX(l,pi,i,nu)*DX(l,k,i,mu)*DFThermalDrOnr(r,t,beta,a);
		}
		else {
			res += (-pow(2.0*PI,2))*DX(l,pk,k,nu)*DX(l,pi,i,nu)*DFThermalDt(r,t,beta,a);
		}			
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PVthrDisjoint_nr
template<uint Dim>
void PVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {	
	ErgVthrDisjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRVthrDisjoint_nr
template<uint Dim>
void PRVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeighDisjoint(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += (-pow(2.0*PI,2))*DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DX(l,k,i,mu)*DFThermalDrOnr(r,t,beta,a);
		}
		else {
			res += (-pow(2.0*PI,2))*DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DFThermalDt(r,t,beta,a);
		}		
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PGaussian_nr
template<uint Dim>
void PGaussian_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v) {
	ErgGaussian_nr(l,j,mu,k,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussian_nr
template<uint Dim>
void PRGaussian_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v) {
	uint pk = posNeigh(k,l.size()), pi;
	number res = 0.0;
	for (uint i=0; i<=j; i++) {
		pi = posNeigh(i,l.size());
		res += (-2.0/a/a)*DX(l,pk,k,nu)*DX(l,pi,i,nu)*DX(l,k,i,mu)*exp(-DistanceSquared(l[i],l[k])/a/a);
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PGaussianDisjoint_nr
template<uint Dim>
void PGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	ErgGaussianDisjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussianDisjoint_nr
template<uint Dim>
void PRGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0;
	for (uint i=0; i<=j; i++) {
		pi = posNeighDisjoint(i,l.size());
		res += (-2.0/a/a)*DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DX(l,k,i,mu)*exp(-DistanceSquared(l[i],l[k])/a/a);
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PGaussianLRDisjoint_nr
template<uint Dim>
void PGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		PGaussianDisjoint_nr(l,j,mu,k,beta,a,f,v);
}

// PRGaussianLRDisjoint_nr
template<uint Dim>
void PRGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		PRGaussianDisjoint_nr(l,j,mu,k,nu,beta,a,f,v);
}

// PGaussianThermal_nr
template<uint Dim>
void PGaussianThermal_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	ErgGaussianThermal_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PGaussianThermal2_nr
template<uint Dim>
void PGaussianThermal2_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	ErgGaussianThermal2_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussianThermal_nr
template<uint Dim>
void PRGaussianThermal_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeigh(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += DX(l,pk,k,nu)*DX(l,pi,i,nu)*DX(l,k,i,mu)*DGThermalDrOnr(r,t,beta,a);
		}
		else {
			res += DX(l,pk,k,nu)*DX(l,pi,i,nu)*DGThermalDt(r,t,beta,a);
		}	
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PRGaussianThermal2_nr
template<uint Dim>
void PRGaussianThermal2_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeigh(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += DX(l,pk,k,nu)*DX(l,pi,i,nu)*DX(l,k,i,mu)*DHThermalDrOnr(r,t,beta,a);
		}
		else {
			res += DX(l,pk,k,nu)*DX(l,pi,i,nu)*DHThermalDt(r,t,beta,a);
		}	
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PGaussianThermalDisjoint_nr
template<uint Dim>
void PGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	ErgGaussianThermalDisjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PGaussianThermal2Disjoint_nr
template<uint Dim>
void PGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	ErgGaussianThermal2Disjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussianThermalDisjoint_nr
template<uint Dim>
void PRGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeighDisjoint(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DX(l,k,i,mu)*DGThermalDrOnr(r,t,beta,a);
		}
		else {
			res += DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DGThermalDt(r,t,beta,a);
		}	
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PRGaussianThermal2Disjoint_nr
template<uint Dim>
void PRGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeighDisjoint(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DX(l,k,i,mu)*DHThermalDrOnr(r,t,beta,a);
		}
		else {
			res += DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DHThermalDt(r,t,beta,a);
		}	
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PGaussianThermalLRDisjoint_nr
template<uint Dim>
void PGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		PGaussianThermalDisjoint_nr(l,j,mu,k,beta,a,f,v);
}

// PGaussianThermal2LRDisjoint_nr
template<uint Dim>
void PGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		PGaussianThermal2Disjoint_nr(l,j,mu,k,beta,a,f,v);
}

// PRGaussianThermalLRDisjoint_nr
template<uint Dim>
void PRGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		PRGaussianThermalDisjoint_nr(l,j,mu,k,nu,beta,a,f,v);
}

// PRGaussianThermal2LRDisjoint_nr
template<uint Dim>
void PRGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		PRGaussianThermal2Disjoint_nr(l,j,mu,k,nu,beta,a,f,v);
}

/*----------------------------------------------------------------------------------------------------------------------------
	2 - loopToVector, vectorToLoop
----------------------------------------------------------------------------------------------------------------------------*/

// loopToVector
template<uint Dim>
void loopToVector(const Loop<Dim>& l, vec& v) {
	v.resize(Dim*l.size());
	uint j,k;
	for (j=0; j<l.size(); j++) {
		for (k=0; k<Dim; k++)
			v[j*Dim+k] = (l[j])[k];
	}
}

// vectorToLoop
template<uint Dim>
void vectorToLoop(const vec& v, Loop<Dim>& l) {
	if (v.size()<Dim*l.size()) {
		cerr << "vectorToLoop error: v.size() = " << v.size() << ", loop needs vector of size = " << Dim*l.size() << endl;
	}
	uint j,k;
	for (j=0; j<l.size(); j++) {
		for (k=0; k<Dim; k++)
			(l[j])[k] = v[j*Dim+k];
	}
}

// printAsLoop
void printAsLoop(const string& f, const uint& Dim, const vec& v) {
	if (v.size()%Dim==0) {
		uint N = v.size()/Dim, mu;
		ofstream os;
		os.open(f.c_str());
		os << left << setprecision(16);
		for (uint j=0; j<N; j++) {
			for (mu=0; mu<Dim; mu++) {
				os << setw(24) << v[j*Dim+mu];
			}
			os << endl;
		}
		os.close();
	}
	else {
		cerr << "printAsLoop error: v.size()%Dim!=0" << endl;
	}
}

// printAsLoop
void printAsLoop(const string& f, const uint& Dim, const vec& v, const uint len) {
	if (len%Dim==0) {
		uint N = len/Dim, mu;
		ofstream os;
		os.open(f.c_str());
		if (!os.good()) {
			cerr << "printAsLoop error: " << f << " not opened properly" << endl;
		}
		os << left << setprecision(16);
		for (uint j=0; j<N; j++) {
			for (mu=0; mu<Dim; mu++) {
				os << setw(24) << v[j*Dim+mu];
			}
			os << endl;
		}
		os.close();
	}
	else {
		cerr << "printAsLoop error: len%Dim!=0" << endl;
	}
}


/*----------------------------------------------------------------------------------------------------------------------------
	3 - filename functions
----------------------------------------------------------------------------------------------------------------------------*/

// filenameLoop
template<uint Dim>
Filename filenameLoopNR(const Parameters& p, const string& base) {
	Filename f = base+"data/nr/loops/dim_"+nts(Dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(p.P4)\
		+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+".dat";
	if (p.Ng>0)
		(f.Extras).push_back(StringPair("Ng",nts(p.Ng)));
	return f;
}

// filenameThermal
template<uint Dim>
Filename filenameThermalNR(const Parameters& p, const string& base) {
	Filename f = base+"data/nr/loops/dim_"+nts(Dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(p.P4)\
		+"_a_"+nts(p.Epsi)+".dat";
	if (abs(p.P4)>MIN_NUMBER)
		(f.Extras).push_back(StringPair("mu",nts(p.Mu)));
	(f.Extras).push_back(StringPair("T",nts(p.T)));
	if (p.Ng>0)
		(f.Extras).push_back(StringPair("Ng",nts(p.Ng)));
	return f;
}

/*----------------------------------------------------------------------------------------------------------------------------
	4 - explicit instantiation
----------------------------------------------------------------------------------------------------------------------------*/

// dim 2
template void L<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void LDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void S0<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void Sm<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void Angle<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void PseudoAngle<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void FGamma<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void InlineCurvatureMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void InlineCurvatureMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void InlineCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void InlineCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void CuspCurvatureMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void CuspCurvatureMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void CuspCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void CuspCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void KGMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void KGMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMaxPlane<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void KGMaxPlane<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMaxDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void KGMaxDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void KGMaxPlaneDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void KGMaxPlaneDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void mdPX_nr<2>(const Loop<2>& l, const uint& loc, const Point<2>& P, const number& f, vec& v);
template void mdL_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void ddL_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void mdS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void ddS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void mdsqrtS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& sqrt4s0, const number& p, vec& v);
template void ddsqrtS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
								 const number& sqrt4s0, const number& p, mat& m);
template void mdVor_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& a, const number& p, vec& v);
template void mdFGamma_nr<2>(const Loop<2>& l, const uint& loc, const number& f, vec& v);
template void ddFGamma_nr<2>(const Loop<2>& l, const uint& loc, const number& f, mat& v);
template void loopToVector<2>(const Loop<2>&,vec&);
template void vectorToLoop<2>(const vec&, Loop<2>&);
template Filename filenameLoopNR<2>(const Parameters& p, const string&);

// mdI0_nr<2>
template <> void mdI0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& f, vec& v) {
	if (mu==1) {
		uint nj = (j==0? (l.size()-1): j-1);
		v[j*2+mu] += -f*((l[nj])[0]-(l[j])[0]);
	}
	else if (mu==0) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		v[j*2+mu] += -f*((l[pj])[1]-(l[j])[1]);
	}
}

// ddI0_nr<2>
template <> void ddI0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<2>& l, const number& f, mat& m) {
	if (mu==1 && nu==0) {
		uint nj = (j==0? (l.size()-1): j-1);
		if (k==nj)
			m(j*2+mu,k*2+nu) += f;
		else if (k==j)
			m(j*2+mu,k*2+nu) -= f;
	}
	else if (mu==0 && nu==1) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		if (k==pj)
			m(j*2+mu,k*2+nu) += f;
		else if (k==j)
			m(j*2+mu,k*2+nu) -= f;
	}
}	

// I0<2>
template <> void I0<2> (const uint& j, const Loop<2>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	result += f*(l[j])[0]*((l[pj])[1]-(l[j])[1]);
}	

// dim 4
template void L<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void LDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void DistPow<4>(const uint& j, const Loop<4>& l, const number& w, const number& f, number& result);
template void S0<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void S0Disjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void I0<4> (const uint& j, const Loop<4>& l, const number& f, number& result);
template void I0Disjoint<4> (const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void In<4> (const uint& j, const Loop<4>& l, const number& n, const number& f, number& result);
template void InDisjoint<4> (const uint& j, const Loop<4>& l, const number& n, const number& beta, const number& f, number& result);
template void Sm<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void Angle<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void FGamma<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void Gaussian<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result);
template void GaussianDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void GaussianLRDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void GaussianThermal<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void GaussianThermalDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void GaussianThermalLRDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void GaussianThermal2<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void GaussianThermal2Disjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void GaussianThermal2LRDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void MaxXn<4>(const uint& j, const uint& k, const Loop<4>& l, const uint& n, const number& f, number& result);
template void MinXnDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const uint& n, const number& f, number& result);
template void InlineCurvatureMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void InlineCurvatureMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void InlineCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void InlineCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void CuspCurvatureMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void CuspCurvatureMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void CuspCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void CuspCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void KGMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void KGMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMaxPlane<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void KGMaxPlane<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMaxDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void KGMaxDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void KGMaxPlaneDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void KGMaxPlaneDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void mdPX_nr<4>(const Loop<4>& l, const uint& loc, const Point<4>& P, const number& f, vec& v);
template void mdL_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void ddL_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void mdLDisjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& p, vec& v);
template void ddLDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
				const number& beta, const number& p, mat& m);
template void mdDistPow_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& w, const number& p, vec& v);
template void ddDistPow_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& w, const number& p, mat& m);
template void mdS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void ddS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void mdS0Disjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& p, vec& v);
template void ddS0Disjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
			const number& p, mat& m);
template void mdsqrtS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& sqrt4s0, const number& p, vec& v);
template void ddsqrtS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
								 const number& sqrt4s0, const number& p, mat& m);
template void mdsqrtS0Disjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& sqrt4s0,\
			 const number& beta, const number& p, vec& v);
template void ddsqrtS0Disjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
								 const number& sqrt4s0, const number& beta, const number& p, mat& m);
template void mdI0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& f, vec& v);
template void ddI0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<4>& l, const number& f, mat& m);	
template void mdI0Disjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& f, vec& v);
template void ddI0Disjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& beta, const number& f, mat& m);
template void mdIn_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& n, const number& f, vec& v);
template void ddIn_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<4>& l, const number& n, const number& f, mat& m);	
template void mdInDisjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& n, \
							const number& beta, const number& f, vec& v);
template void ddInDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& n,\
 							const number& beta, const number& f, mat& m);
template void mdVor_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& a, const number& p, vec& v);
template void mdFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, vec& v);
template void ddFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, mat& v);
template void ErgS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void ErgS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void ErgL_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void ErgLDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void ErgsqrtS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, number& erg);
template void ErgsqrtS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, number& erg);
template void ErgI0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void ErgI0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void ErgIn_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& f, number& erg);
template void ErgInDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, number& erg);
template void ErgVor_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg);
template void ErgVthr_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgVthrDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgVthrDisjointLR_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussian_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, number& erg);
template void ErgGaussianDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussianLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussianThermal_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussianThermalDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussianThermalLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussianThermal2_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussianThermal2Disjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void ErgGaussianThermal2LRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void PS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void PS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void PL_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void PLDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void PsqrtS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, vec& v);
template void PsqrtS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, vec& v);
template void PI0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void PI0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void PIn_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& f, vec& v);
template void PInDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, vec& v);
template void PVor_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v);
template void PRVor_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);
template void PVthr_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRVthr_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PVthrDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRVthrDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussian_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, vec& v);
template void PRGaussian_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);
template void PGaussianDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussianLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussianThermal_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianThermal_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussianThermalDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianThermalDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussianThermalLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianThermalLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussianThermal2_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianThermal2_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussianThermal2Disjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianThermal2Disjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void PGaussianThermal2LRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void PRGaussianThermal2LRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void loopToVector<4>(const Loop<4>&,vec&);
template void vectorToLoop<4>(const vec&, Loop<4>&);
template Filename filenameLoopNR<4>(const Parameters& p, const string&);
template Filename filenameThermalNR<4>(const Parameters& p, const string&);

template <> void PseudoAngle<4>(const uint& j, const Loop<4>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	number sin_gamma;
	
	sin_gamma = (l[pj]-l[j])[2]*(l[j]-l[nj])[3] - (l[pj]-l[j])[3]*(l[j]-l[nj])[2];
	sin_gamma /= sqrt((pow((l[pj]-l[j])[2],2)+pow((l[pj]-l[j])[3],2))*(pow((l[j]-l[nj])[2],2)+pow((l[j]-l[nj])[3],2)));
	
	result += f*(asin(sin_gamma));
}

// Vor
template <> void Vor<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])/(DistanceSquared(l[j],l[k])+a*a);	
	}
}

// Vlr
template <> void Vlr<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])/(0.25*DistanceSquared(l[pj]+l[j],l[pk]+l[k])+a*a);	
	}
}

// Ver
template <> void Ver<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
		number d2 = DistanceSquared(l[j],l[k]);
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])*exp(-a*a/d2)/d2;	
	}
}

// Vdr
template <> void Vdr<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])/pow(DistanceSquared(l[j],l[k]),1.0-a/2.0);	
	}
}

// Vthr
template <> void Vthr<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeigh(j,l.size());
		uint pk = posNeigh(k,l.size());
		
		number r = SpatialDistance(l[j],l[k]);
		number t = DX(l,k,j,3);
	
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*FThermal(r,t,beta,a);	
	}
}

// VthrDisjoint
template <>
void VthrDisjoint<4> (const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result) {
	if (k<=j) {
		uint pj = posNeighDisjoint(j,l.size());
		uint pk = posNeighDisjoint(k,l.size());
		
		number r = SpatialDistance(l[j],l[k]);
		number t = DXDisjoint(l,k,j,3,beta);
	
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*DotDisjoint(l[pj],l[j],l[pk],l[k],beta)*FThermal(r,t,beta,a);	
	}
}

// VthrDisjointLR
template <>
void VthrDisjointLR<4> (const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		VthrDisjoint(j,k,l,beta,a,f,result);
	else
		VthrDisjoint(j,k,l,beta,0.0,f,result);
}

// VnonrelDisjoint
template <>
void VnonrelDisjoint<4> (const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result) {
	if (j<l.size()/2) {
		number dt = 2.0*beta/(number)l.size();
		uint oj = oppNeigh(j,l.size());

		result += f*dt/SpatialDistance(l[j],l[oj]); // dt/r, not lorentz invariant but galilean invariant
	}
}

// VnonrelrDisjoint
template <>
void VnonrelrDisjoint<4> (const uint& j, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result) {
	if (j<l.size()/2) {
		uint oj = oppNeigh(j,l.size());
		number dt = 2.0*beta/(number)l.size();
		number r2 = SpatialDistanceSquared(l[j],l[oj]);
		number a2 = pow(a,2);

		result += f*dt*(1.0/sqrt(r2 + a2) - exp(-r2/a2)/a); // not lorentz invariant but galilean invariant
	}
}

// mdVor_nr
template <> void mdVor_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v) {
	number res = 0.0;
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	uint pi = (i==(l.size()-1)? 0: i+1);
		
	if (i!=j) {
		number B_ij = a*a + DistanceSquared(l[i],l[j]);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		res += 2.0*DX(l,i,pi,mu)/B_ij - 4.0*DX(l,j,i,mu)*T_ij/pow(B_ij,2);
	}
	
	if (i!=mj) {
		number B_imj = a*a + DistanceSquared(l[i],l[mj]);
		res += -2.0*DX(l,i,pi,mu)/B_imj;
	}
	
	//coincident terms
	if (i==j)
		res += 2.0*(l[j])[mu]/a/a;
	if (i==mj)
		res += -(l[mj])[mu]/a/a;
	if (i==pj)
		res += -(l[pj])[mu]/a/a;

	v[j*4+mu] += -f*res;
}

// mdVlr_nr
template <> void mdVlr_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v) {
	number B_ij, B_imj, T_ij, T_imj, res = 0.0;
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	uint pi = (i==(l.size()-1)? 0: i+1);
		
	if (i!=j) {
		B_ij = a*a + 0.25*DistanceSquared(l[pi]+l[i],l[pj]+l[j]);
		T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		res =  + 2.0*DX(l,i,pi,mu)/B_ij \
				+ ( (-DX(l,j,i,mu) - DX(l,pj,pi,mu)) \
				- (DX(l,j,i,mu) + DX(l,pj,pi,mu)) )*T_ij/(2.0*pow(B_ij,2));
	}
	if (i!=mj) {	
		B_imj = a*a + 0.25*DistanceSquared(l[pi]+l[i],l[j]+l[mj]);
		T_imj = Dot(l[pi],l[i],l[j],l[mj]);
		res +=  -2.0*DX(l,i,pi,mu)/B_imj 
				+( - (-DX(l,mj,i,mu) - DX(l,j,pi,mu)) \
				+ (DX(l,mj,i,mu) + DX(l,j,pi,mu)) )*T_imj/(2.0*pow(B_imj,2));
	}
		
	//coincident terms
	if (i==j)
		res += 2.0*(l[j])[mu]/a/a;
	if (i==mj)
		res += -(l[mj])[mu]/a/a;
	if (i==pj)
		res += -(l[pj])[mu]/a/a;

	v[j*4+mu] += -f*res;
}

// mdVer_nr
template <> void mdVer_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v) {
	number res = 0.0;
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	uint pi = (i==(l.size()-1)? 0: i+1);
		
	if (i!=j) {
		number B_ij = DistanceSquared(l[i],l[j]);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		number E_ij = exp(-a*a/B_ij);
		res += (+ (2.0*DX(l,i,pi,mu))/B_ij \
 				+ (4.0*pow(a,2)*DX(l,j,i,mu)*T_ij)/pow(B_ij,3) \
 				- (4.0*DX(l,j,i,mu)*T_ij)/pow(B_ij,2) )*E_ij;
	}
	if (i!=mj) {
		number B_imj = DistanceSquared(l[i],l[mj]);
		number E_imj = exp(-a*a/B_imj);
		res += (-2.0*DX(l,i,pi,mu)/B_imj)*E_imj;
	}

	v[j*4+mu] += -f*res;
}

// mdVdr_nr
template <> void mdVdr_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v) {
	number res = 0.0;
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	uint pi = (i==(l.size()-1)? 0: i+1);
		
	if (i!=j) {
		number B_ij = DistanceSquared(l[i],l[j]);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		res +=   2.0*pow(B_ij,-1.0 + a/2.0)*DX(l,i,pi,mu) \
 				+  2.0*(-2.0 + a)*pow(B_ij,-2.0 + a/2.0)*DX(l,j,i,mu)*T_ij;
	}
	if (i!=mj) {
		number B_imj = DistanceSquared(l[i],l[mj]);
		res += - 2.0*(pow(B_imj,-1.0 + a/2.0)*DX(l,i,pi,mu));
	}

	v[j*4+mu] += -f*res;
}

// mdVthr_nr
template <> void mdVthr_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
			const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DX(l,j,i,3); // checked order
		number FThermal_ij = FThermal(r_ij,t_ij,beta,a);
		number DFThermalDrOnr_ij = DFThermalDrOnr(r_ij,t_ij,beta,a);
		number DFThermalDt_ij = DFThermalDt(r_ij,t_ij,beta,a);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		
		res += 2.0*FThermal_ij*(-DX(l,i,mu));
		if (mu<3)
			res += 2.0*DFThermalDrOnr_ij*DX(l,j,i,mu)*T_ij;
		else
			res += 2.0*DFThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DX(l,mj,i,3); // checked order
		number FThermal_imj = FThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*FThermal_imj*(-DX(l,i,mu)); //
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (i==j)
		res += (-1.0/pow(2.0*PI,2))*2.0*(l[j])[mu]/a/a; 
	if (i==mj)
		res += (-1.0/pow(2.0*PI,2))*(-(l[mj])[mu]/a/a);
	if (i==pj)
		res += (-1.0/pow(2.0*PI,2))*(-(l[pj])[mu]/a/a);
		
	v[j*4+mu] += -f*(-pow(2.0*PI,2))*res;
}

// mdVthrDisjoint_nr
template <> void mdVthrDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
			const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeighDisjoint(j,l.size());
	uint mj = negNeighDisjoint(j,l.size());
	uint pi = posNeighDisjoint(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DXDisjoint(l,j,i,3,beta); // checked order
		number FThermal_ij = FThermal(r_ij,t_ij,beta,a);
		number DFThermalDrOnr_ij = DFThermalDrOnr(r_ij,t_ij,beta,a);
		number DFThermalDt_ij = DFThermalDt(r_ij,t_ij,beta,a);
		number T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
		
		res += 2.0*FThermal_ij*(-DXDisjoint(l,i,mu,beta));
		if (mu<3)
			res += 2.0*DFThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij;
		else
			res += 2.0*DFThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DXDisjoint(l,mj,i,3,beta); // checked order
		number FThermal_imj = FThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*FThermal_imj*(-DXDisjoint(l,i,mu,beta)); //
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (mu==3) {
		if (i==j) { // bit of a fudge in terms of sum over i
			res += (-1.0/pow(2.0*PI,2))*(DXDisjoint(l,j,mj,mu,beta) + DXDisjoint(l,j,pj,mu,beta))/a/a;
		}
	}
	else {
		if (i==j)
			res += (-1.0/pow(2.0*PI,2))*2.0*(l[j])[mu]/a/a;
		if (i==mj)
			res += (-1.0/pow(2.0*PI,2))*(-(l[mj])[mu]/a/a);
		if (i==pj)
			res += (-1.0/pow(2.0*PI,2))*(-(l[pj])[mu]/a/a);
	}
		
	v[j*4+mu] += -f*(-pow(2.0*PI,2))*res;
}

// mdVthrDisjointLR_nr
template <> void mdVthrDisjointLR_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
			const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && i<N/2) || (j>=N/2 && i>=N/2))
		mdVthrDisjoint_nr(j,mu,i,l,beta,a,f,v);
	else
		mdVthrDisjoint_nr(j,mu,i,l,beta,0.0,f,v);
}


// mdVnonrelDisjoint_nr
template <>
void mdVnonrelDisjoint_nr<4> (const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& f, vec& v) {
	if (mu<3) {
		uint oj = oppNeigh(j,l.size());
		number dt = 2.0*beta/(number)l.size();
		v[j*4+mu] += -f*dt*DX(l,j,oj,mu)*(-1.0/pow(SpatialDistance(l[j],l[oj]),3));
	}
}

// mdVnonrelrDisjoint_nr
template <>
void mdVnonrelrDisjoint_nr<4> (const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& a, const number& f, vec& v) {
	if (mu<3) {
		uint oj = oppNeigh(j,l.size());
		number dt = 2.0*beta/(number)l.size();
		number r2 = SpatialDistanceSquared(l[j],l[oj]);
		number a2 = pow(a,2);
		v[j*4+mu] += -f*dt*DX(l,j,oj,mu)*(-1.0/pow(r2 + a2,3.0/2.0) + 2.0*exp(-r2/a2)/a/a2 );
	}
}

// mdGaussian_nr
template <> void mdGaussian_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v) {
 
	number res = 0.0;
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number B_ij = DistanceSquared(l[i],l[j]);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		number E_ij = exp(-B_ij/a/a);
		res += + 2.0*E_ij*(-DX(l,i,mu)) \
 			- (4.0*E_ij*DX(l,j,i,mu)*T_ij)/pow(a,2);
	}
	if (i!=mj) {
		number B_imj = DistanceSquared(l[i],l[mj]);
		number E_imj = exp(-B_imj/a/a);
		res += - 2.0*E_imj*(-DX(l,i,mu));
	}
	
	//coincident terms
	if (i==j)
		res += 2.0*(l[j])[mu];
	if (i==mj)
		res += -(l[mj])[mu];
	if (i==pj)
		res += -(l[pj])[mu];

	v[j*4+mu] += -f*res;

}

// mdGaussianDisjoint_nr
template <> void mdGaussianDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v) {
				
 	number res = 0.0;
	uint pj = posNeighDisjoint(j,l.size());
	uint mj = negNeighDisjoint(j,l.size());
	uint pi = posNeighDisjoint(i,l.size());
		
	if (i!=j) {
		number B_ij = DistanceSquaredDisjoint(l[i],l[j],beta);
		number T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
		number E_ij = exp(-B_ij/a/a);
		res += + 2.0*E_ij*(-DXDisjoint(l,i,mu,beta)) \
 			- (4.0*E_ij*DXDisjoint(l,j,i,mu,beta)*T_ij)/pow(a,2);
	}
	if (i!=mj) {
		number B_imj = DistanceSquaredDisjoint(l[i],l[mj],beta);
		number E_imj = exp(-B_imj/a/a);
		res += - 2.0*E_imj*(-DXDisjoint(l,i,mu,beta));
	}
		
	//coincident terms	
	if (mu==3) {
		if (i==j) { // bit of a fudge in terms of sum over i
			res += (DXDisjoint(l,j,mj,mu,beta) + DXDisjoint(l,j,pj,mu,beta));
		}
	}
	else {
		if (i==j)
			res += 2.0*(l[j])[mu];
		if (i==mj)
			res += -(l[mj])[mu];
		if (i==pj)
			res += -(l[pj])[mu];
	}

	v[j*4+mu] += -f*res;	

}

// mdGaussianLRDisjoint_nr
template <> void mdGaussianLRDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && i<N/2) || (j>=N/2 && i>=N/2))
		mdGaussianDisjoint_nr(j,mu,i,l,beta,a,f,v);
}

// mdGaussianThermal_nr
template <> void mdGaussianThermal_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
				const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DX(l,j,i,3); // checked order
		number GThermal_ij = GThermal(r_ij,t_ij,beta,a);
		number DGThermalDrOnr_ij = DGThermalDrOnr(r_ij,t_ij,beta,a);
		number DGThermalDt_ij = DGThermalDt(r_ij,t_ij,beta,a);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		
		res += 2.0*GThermal_ij*(-DX(l,i,mu));
		if (mu<3)
			res += 2.0*DGThermalDrOnr_ij*DX(l,j,i,mu)*T_ij;
		else
			res += 2.0*DGThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DX(l,mj,i,3); // checked order
		number GThermal_imj = GThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*GThermal_imj*(-DX(l,i,mu)); //
	}
	
	//coincident terms
	if (i==j)
		res += 2.0*(l[j])[mu]; 
	if (i==mj)
		res += -(l[mj])[mu];
	if (i==pj)
		res += -(l[pj])[mu];
		
	v[j*4+mu] += -f*res;
}

// mdGaussianThermal2_nr
template <> void mdGaussianThermal2_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
				const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DX(l,j,i,3); // checked order
		number HThermal_ij = HThermal(r_ij,t_ij,beta,a);
		number DHThermalDrOnr_ij = DHThermalDrOnr(r_ij,t_ij,beta,a);
		number DHThermalDt_ij = DHThermalDt(r_ij,t_ij,beta,a);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		
		res += 2.0*HThermal_ij*(-DX(l,i,mu));
		if (mu<3)
			res += 2.0*DHThermalDrOnr_ij*DX(l,j,i,mu)*T_ij;
		else
			res += 2.0*DHThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DX(l,mj,i,3); // checked order
		number HThermal_imj = HThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*HThermal_imj*(-DX(l,i,mu)); //
	}
	
	//coincident terms
	if (i==j)
		res += 2.0*(l[j])[mu]; 
	if (i==mj)
		res += -(l[mj])[mu];
	if (i==pj)
		res += -(l[pj])[mu];
		
	v[j*4+mu] += -f*res;
}

// mdGaussianThermalDisjoint_nr
template <> void mdGaussianThermalDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v) {			
 	number res = 0.0;
	
	uint pj = posNeighDisjoint(j,l.size());
	uint mj = negNeighDisjoint(j,l.size());
	uint pi = posNeighDisjoint(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DXDisjoint(l,j,i,3,beta); // checked order
		number GThermal_ij = GThermal(r_ij,t_ij,beta,a);
		number DGThermalDrOnr_ij = DGThermalDrOnr(r_ij,t_ij,beta,a);
		number DGThermalDt_ij = DGThermalDt(r_ij,t_ij,beta,a);
		number T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
		
		res += 2.0*GThermal_ij*(-DXDisjoint(l,i,mu,beta));
		if (mu<3)
			res += 2.0*DGThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij;
		else
			res += 2.0*DGThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DXDisjoint(l,mj,i,3,beta); // checked order
		number GThermal_imj = GThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*GThermal_imj*(-DXDisjoint(l,i,mu,beta)); //
	}
	
	//coincident terms	
	if (mu==3) {
		if (i==j) { // bit of a fudge in terms of sum over i
			res += (DXDisjoint(l,j,mj,mu,beta) + DXDisjoint(l,j,pj,mu,beta));
		}
	}
	else {
		if (i==j)
			res += 2.0*(l[j])[mu];
		if (i==mj)
			res += -(l[mj])[mu];
		if (i==pj)
			res += -(l[pj])[mu];
	}
		
	v[j*4+mu] += -f*res;
}


// mdGaussianThermal2Disjoint_nr
template <> void mdGaussianThermal2Disjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v) {			
 	number res = 0.0;
	
	uint pj = posNeighDisjoint(j,l.size());
	uint mj = negNeighDisjoint(j,l.size());
	uint pi = posNeighDisjoint(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DXDisjoint(l,j,i,3,beta); // checked order
		number HThermal_ij = HThermal(r_ij,t_ij,beta,a);
		number DHThermalDrOnr_ij = DHThermalDrOnr(r_ij,t_ij,beta,a);
		number DHThermalDt_ij = DHThermalDt(r_ij,t_ij,beta,a);
		number T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
		
		res += 2.0*HThermal_ij*(-DXDisjoint(l,i,mu,beta));
		if (mu<3)
			res += 2.0*DHThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij;
		else
			res += 2.0*DHThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DXDisjoint(l,mj,i,3,beta); // checked order
		number HThermal_imj = HThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*HThermal_imj*(-DXDisjoint(l,i,mu,beta)); //
	}
	
	//coincident terms	
	if (mu==3) {
		if (i==j) { // bit of a fudge in terms of sum over i
			res += (DXDisjoint(l,j,mj,mu,beta) + DXDisjoint(l,j,pj,mu,beta));
		}
	}
	else {
		if (i==j)
			res += 2.0*(l[j])[mu];
		if (i==mj)
			res += -(l[mj])[mu];
		if (i==pj)
			res += -(l[pj])[mu];
	}
		
	v[j*4+mu] += -f*res;
}

// mdGaussianThermalLRDisjoint_nr
template <> void mdGaussianThermalLRDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v) {		
	uint N = l.size();
	if ( (j<N/2 && i<N/2) || (j>=N/2 && i>=N/2) )
			mdGaussianThermalDisjoint_nr(j,mu,i,l,beta,a,f,v);
}

// mdGaussianThermal2LRDisjoint_nr
template <> void mdGaussianThermal2LRDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v) {		
	uint N = l.size();
	if ( (j<N/2 && i<N/2) || (j>=N/2 && i>=N/2) )
			mdGaussianThermal2Disjoint_nr(j,mu,i,l,beta,a,f,v);
}

// ddVor_nr
template <> void ddVor_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = (j==0? (l.size()-1): j-1);
	uint pj = (j==(l.size()-1)? 0: j+1);		
	uint mk = (k==0? (l.size()-1): k-1);
	uint pk = (k==(l.size()-1)? 0: k+1);
	
	number B_jk = a*a + DistanceSquared(l[j],l[k]);
	number B_mjk = a*a + DistanceSquared(l[mj],l[k]);
	number B_jmk = a*a + DistanceSquared(l[j],l[mk]);
	number B_mjmk = a*a + DistanceSquared(l[mj],l[mk]);
	
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 	 + 2.0/B_jk \
					 + 2.0/B_mjmk \
					 + 4*T_jk/pow(B_jk,2);
		if (k!=mj)
			res += - 2.0/B_mjk; //
		if (k!=pj)
			res += - 2.0/B_jmk; //
	}
	
	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res +=  - 4*DX(l,j,pj,nu)*DX(l,j,k,mu)/pow(B_jk,2) \
				 + 4*DX(l,j,k,nu)*DX(l,k,pk,mu)/pow(B_jk,2) \
				 - 16*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk/pow(B_jk,3); //
	if (k!=mj)
		res += - 4*DX(l,mj,k,nu)*DX(l,k,pk,mu)/pow(B_mjk,2); //
	if (k!=pj)
		res += + 4*DX(l,j,pj,nu)*DX(l,j,mk,mu)/pow(B_jmk,2); //
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number B_ij, B_imj, T_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = (i==(l.size()-1)? 0: i+1);
			B_ij = a*a + DistanceSquared(l[i],l[j]);
			B_imj = a*a + DistanceSquared(l[i],l[mj]);
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);
			
			if (k==j && i!=j) {
				res += 	 - 4.0*DX(l,i,pi,nu)*DX(l,j,i,mu)/pow(B_ij,2) \
						 - 4.0*DX(l,i,pi,mu)*DX(l,j,i,nu)/pow(B_ij,2) \
						 + 16.0*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij/pow(B_ij,3);
				if (mu==nu)
					res += - 4.0*T_ij/pow(B_ij,2); //
			}
			if (k==mj && i!=mj) 
				res += 4.0*DX(l,i,pi,mu)*DX(l,mj,i,nu)/pow(B_imj,2);
			if (k==pj && i!=j) 
				res += 4.0*DX(l,i,pi,nu)*DX(l,j,i,mu)/pow(B_ij,2);
			
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0/a/a;
	if (k==mj && mu==nu)
		res += -1.0/a/a;
	if (k==pj && mu==nu)
		res += -1.0/a/a;
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// ddVlr_nr
template <> void ddVlr_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = (j==0? (l.size()-1): j-1);
	uint pj = (j==(l.size()-1)? 0: j+1);		
	uint mk = (k==0? (l.size()-1): k-1);
	uint pk = (k==(l.size()-1)? 0: k+1);
	
	number B_jk = a*a + 0.25*DistanceSquared(l[pj]+l[j],l[pk]+l[k]);
	number B_mjk = a*a + 0.25*DistanceSquared(l[j]+l[mj],l[pk]+l[k]);
	number B_jmk = a*a + 0.25*DistanceSquared(l[pj]+l[j],l[k]+l[mk]);
	number B_mjmk = a*a + 0.25*DistanceSquared(l[j]+l[mj],l[k]+l[mk]);
	
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	number T_mjk = Dot(l[j],l[mj],l[pk],l[k]);
	number T_jmk = Dot(l[pj],l[j],l[k],l[mk]);
	number T_mjmk = Dot(l[j],l[mj],l[k],l[mk]);
	
	// terms where mu==nu, without sums		
	if (mu==nu) {
		if (k!=j) {
			res += 	+ 2.0/B_mjmk \
					+ 2.0/B_jk \
					+ T_mjmk/pow(B_mjmk,2) \
					+ T_jk/pow(B_jk,2);
		}
		if (k!=mj)
			res += 	- 2.0/B_mjk \
					- T_mjk/pow(B_mjk,2);
		if (k!=pj)
			res += 	- 2.0/B_jmk \
					- T_jmk/pow(B_jmk,2);
	}
	
// terms where mu not nexcessarily equal to nu, without sums#
	if (k!=j) {
		res += - (DX(l,j,pj,nu)*DX(l,j,k,mu))/pow(B_jk,2) \
				 - (DX(l,j,pj,nu)*DX(l,pj,pk,mu))/pow(B_jk,2) \
				 + (DX(l,j,k,nu)*DX(l,k,pk,mu))/pow(B_jk,2) \
				 + (DX(l,pj,pk,nu)*DX(l,k,pk,mu))/pow(B_jk,2) \
				 - (DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(B_jk,3) \
				 - (DX(l,j,k,nu)*DX(l,pj,pk,mu)*T_jk)/pow(B_jk,3) \
				 - (DX(l,j,k,mu)*DX(l,pj,pk,nu)*T_jk)/pow(B_jk,3) \
				 - (DX(l,pj,pk,mu)*DX(l,pj,pk,nu)*T_jk)/pow(B_jk,3) \
				  - (DX(l,mj,mk,mu)*DX(l,j,mj,nu))/pow(B_mjmk,2) \
				 - (DX(l,j,mj,nu)*DX(l,j,k,mu))/pow(B_mjmk,2) \
				 + (DX(l,mj,mk,nu)*DX(l,k,mk,mu))/pow(B_mjmk,2) \
				 + (DX(l,j,k,nu)*DX(l,k,mk,mu))/pow(B_mjmk,2) \
				 - (DX(l,mj,mk,mu)*DX(l,mj,mk,nu)*T_mjmk)/pow(B_mjmk,3) \
				 - (DX(l,mj,mk,nu)*DX(l,j,k,mu)*T_mjmk)/pow(B_mjmk,3) \
				 - (DX(l,mj,mk,mu)*DX(l,j,k,nu)*T_mjmk)/pow(B_mjmk,3) \
				 - (DX(l,j,k,mu)*DX(l,j,k,nu)*T_mjmk)/pow(B_mjmk,3);
	}
	if (k!=mj) {
		res += 	+ (DX(l,mj,k,mu)*DX(l,j,mj,nu))/pow(B_mjk,2) \
				 + (DX(l,j,mj,nu)*DX(l,j,pk,mu))/pow(B_mjk,2) \
				 - (DX(l,mj,k,nu)*DX(l,k,pk,mu))/pow(B_mjk,2) \
				 - (DX(l,j,pk,nu)*DX(l,k,pk,mu))/pow(B_mjk,2) \
				 + (DX(l,mj,k,mu)*DX(l,mj,k,nu)*T_mjk)/pow(B_mjk,3) \
				 + (DX(l,mj,k,nu)*DX(l,j,pk,mu)*T_mjk)/pow(B_mjk,3) \
				 + (DX(l,mj,k,mu)*DX(l,j,pk,nu)*T_mjk)/pow(B_mjk,3) \
				 + (DX(l,j,pk,mu)*DX(l,j,pk,nu)*T_mjk)/pow(B_mjk,3);
	}
	if (k!=pj) {
		res +=    + (DX(l,j,pj,nu)*DX(l,j,mk,mu))/pow(B_jmk,2) \
				+ (DX(l,j,pj,nu)*DX(l,pj,k,mu))/pow(B_jmk,2) \
				- (DX(l,j,mk,nu)*DX(l,k,mk,mu))/pow(B_jmk,2) \
				- (DX(l,pj,k,nu)*DX(l,k,mk,mu))/pow(B_jmk,2) \
				+ (DX(l,j,mk,mu)*DX(l,j,mk,nu)*T_jmk)/pow(B_jmk,3) \
				+ (DX(l,j,mk,nu)*DX(l,pj,k,mu)*T_jmk)/pow(B_jmk,3) \
				+ (DX(l,j,mk,mu)*DX(l,pj,k,nu)*T_jmk)/pow(B_jmk,3) \
				+ (DX(l,pj,k,mu)*DX(l,pj,k,nu)*T_jmk)/pow(B_jmk,3);
	}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number B_ij, B_imj, T_ij, T_imj;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = (i==(l.size()-1)? 0: i+1);		
			
			B_ij = a*a + 0.25*DistanceSquared(l[pj]+l[j],l[pi]+l[i]);
			B_imj = a*a + 0.25*DistanceSquared(l[j]+l[mj],l[pi]+l[i]);
			T_ij = Dot(l[pj],l[j],l[pi],l[i]);
			T_imj = Dot(l[j],l[mj],l[pi],l[i]);
			
			if (k==j && i!=j) {
				res += - (DX(l,i,pi,nu)*(-DX(l,mj,i,mu) - DX(l,j,pi,mu)))/(2.0*pow(B_imj,2)) \
						+ (DX(l,i,pi,nu)*(DX(l,mj,i,mu) + DX(l,j,pi,mu)))/(2.0*pow(B_imj,2)) \
						- (DX(l,i,pi,mu)*(-DX(l,mj,i,nu) - DX(l,j,pi,nu)))/pow(B_imj,2) \
						+ (DX(l,i,pi,nu)*(-DX(l,j,i,mu) - DX(l,pj,pi,mu)))/(2.0*pow(B_ij,2)) \
						- (DX(l,i,pi,nu)*(DX(l,j,i,mu) + DX(l,pj,pi,mu)))/(2.0*pow(B_ij,2)) \
						+ (DX(l,i,pi,mu)*(-DX(l,j,i,nu) - DX(l,pj,pi,nu)))/pow(B_ij,2) \
						- ((-DX(l,mj,i,mu) - DX(l,j,pi,mu))*(-DX(l,mj,i,nu) - DX(l,j,pi,nu))*T_imj)/ \
						(2.0*pow(B_imj,3)) + ((DX(l,mj,i,mu) + DX(l,j,pi,mu))*(-DX(l,mj,i,nu) - DX(l,j,pi,nu))* \
						T_imj)/(2.0*pow(B_imj,3)) \
						+ ((-DX(l,j,i,mu) - DX(l,pj,pi,mu))*(-DX(l,j,i,nu) - DX(l,pj,pi,nu))*T_ij)/ \
						(2.0*pow(B_ij,3)) - ((DX(l,j,i,mu) + DX(l,pj,pi,mu))*(-DX(l,j,i,nu) - DX(l,pj,pi,nu))* \
						T_ij)/(2.0*pow(B_ij,3));		

				if (mu==nu)
					res += + T_imj/pow(B_imj,2) \
							- T_ij/pow(B_ij,2);
			}
			if (k==mj && i!=mj) {
				res +=  (DX(l,i,pi,nu)*(-DX(l,mj,i,mu) - DX(l,j,pi,mu)))/(2.0*pow(B_imj,2)) \
						- (DX(l,i,pi,nu)*(DX(l,mj,i,mu) + DX(l,j,pi,mu)))/(2.0*pow(B_imj,2)) \
						- (DX(l,i,pi,mu)*(-DX(l,mj,i,nu) - DX(l,j,pi,nu)))/pow(B_imj,2) \
						- ((-DX(l,mj,i,mu) - DX(l,j,pi,mu))*(-DX(l,mj,i,nu) - DX(l,j,pi,nu))*T_imj)/ \
						(2.0*pow(B_imj,3)) + ((DX(l,mj,i,mu) + DX(l,j,pi,mu))*(-DX(l,mj,i,nu) - DX(l,j,pi,nu))* \
						T_imj)/(2.0*pow(B_imj,3));
				if (mu==nu)
					res += + T_imj/pow(B_imj,2);
			}
			if (k==pj && i!=j) {
				res += 	 - (DX(l,i,pi,nu)*(-DX(l,j,i,mu) - DX(l,pj,pi,mu)))/(2.0*pow(B_ij,2)) \
						+ (DX(l,i,pi,nu)*(DX(l,j,i,mu) + DX(l,pj,pi,mu)))/(2.0*pow(B_ij,2)) \
						+ (DX(l,i,pi,mu)*(-DX(l,j,i,nu) - DX(l,pj,pi,nu)))/pow(B_ij,2) \
						+ ((-DX(l,j,i,mu) - DX(l,pj,pi,mu))*(-DX(l,j,i,nu) - DX(l,pj,pi,nu))*T_ij)/ \
						(2.0*pow(B_ij,3)) - ((DX(l,j,i,mu) + DX(l,pj,pi,mu))*(-DX(l,j,i,nu) - DX(l,pj,pi,nu))* \
						T_ij)/(2.0*pow(B_ij,3));
				if (mu==nu)
					res += - T_ij/pow(B_ij,2);
			}
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0/a/a;
	if (k==mj && mu==nu)
		res += -1.0/a/a;
	if (k==pj && mu==nu)
		res += -1.0/a/a;
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// ddVer_nr
template <> void ddVer_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = (j==0? (l.size()-1): j-1);
	uint pj = (j==(l.size()-1)? 0: j+1);		
	uint mk = (k==0? (l.size()-1): k-1);
	uint pk = (k==(l.size()-1)? 0: k+1);
	
	number B_jk = DistanceSquared(l[j],l[k]);
	number B_mjk = DistanceSquared(l[mj],l[k]);
	number B_jmk = DistanceSquared(l[j],l[mk]);
	number B_mjmk = DistanceSquared(l[mj],l[mk]);
	
	number E_jk = exp(-a*a/B_jk);
	number E_mjk = exp(-a*a/B_mjk);
	number E_jmk = exp(-a*a/B_jmk);
	number E_mjmk = exp(-a*a/B_mjmk);
	
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 	 + (2.0*E_mjmk)/B_mjmk \
					 + (2.0*E_jk)/B_jk \
					 - (4.0*pow(a,2)*E_jk*T_jk)/pow(B_jk,3) \
					 + (4.0*E_jk*T_jk)/pow(B_jk,2);
		if (k!=mj)
			res +=  - (2.0*E_mjk)/B_mjk; //
		if (k!=pj)
			res += - (2.0*E_jmk)/B_jmk; //
	}
	
	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res +=   + (4.0*pow(a,2)*E_jk*DX(l,j,pj,nu)*DX(l,j,k,mu))/pow(B_jk,3) \
 				- (4.0*E_jk*DX(l,j,pj,nu)*DX(l,j,k,mu))/pow(B_jk,2) \
 				- (4.0*pow(a,2)*E_jk*DX(l,j,k,nu)*DX(l,k,pk,mu))/pow(B_jk,3) \
				 + (4.0*E_jk*DX(l,j,k,nu)*DX(l,k,pk,mu))/pow(B_jk,2) \
				 - (8.0*pow(a,4)*E_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(B_jk,5) \
				 + (32.0*pow(a,2)*E_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(B_jk,4) \
				 - (16.0*E_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(B_jk,3); //
	if (k!=mj)
		res +=  + (4.0*pow(a,2)*E_mjk*DX(l,mj,k,nu)*DX(l,k,pk,mu))/pow(B_mjk,3) \
 				- (4.0*E_mjk*DX(l,mj,k,nu)*DX(l,k,pk,mu))/pow(B_mjk,2); //
	if (k!=pj)
		res += - (4.0*pow(a,2)*E_jmk*DX(l,j,pj,nu)*DX(l,j,mk,mu))/pow(B_jmk,3) \
 				+ (4.0*E_jmk*DX(l,j,pj,nu)*DX(l,j,mk,mu))/pow(B_jmk,2) ; //
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number B_ij, B_imj, E_ij, E_imj, T_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = (i==(l.size()-1)? 0: i+1);
			B_ij = DistanceSquared(l[i],l[j]);
			B_imj = DistanceSquared(l[i],l[mj]);
			E_ij = exp(-a*a/B_ij);
			E_imj = exp(-a*a/B_imj);
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);
			
			if (k==j && i!=j) {
				res += 	 + (4.0*pow(a,2)*E_ij*DX(l,i,pi,nu)*DX(l,j,i,mu))/pow(B_ij,3) \
						 - (4.0*E_ij*DX(l,i,pi,nu)*DX(l,j,i,mu))/pow(B_ij,2) \
						 + (4.0*pow(a,2)*E_ij*DX(l,i,pi,mu)*DX(l,j,i,nu))/pow(B_ij,3) \
						 - (4.0*E_ij*DX(l,i,pi,mu)*DX(l,j,i,nu))/pow(B_ij,2) \
						 + (8.0*pow(a,4)*E_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(B_ij,5) \
						 - (32.0*pow(a,2)*E_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(B_ij,4) \
						 + (16.0*E_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(B_ij,3);
				if (mu==nu)
					res +=  + (4.0*pow(a,2)*E_ij*T_ij)/pow(B_ij,3) \
 							- (4.0*E_ij*T_ij)/pow(B_ij,2); //
			}
			if (k==mj && i!=mj) 
				res +=  (-4.0*pow(a,2)*E_imj*DX(l,i,pi,mu)*DX(l,mj,i,nu))/pow(B_imj,3) \
 						+ (4.0*E_imj*DX(l,i,pi,mu)*DX(l,mj,i,nu))/pow(B_imj,2);
			if (k==pj && i!=j) 
				res += - (4.0*pow(a,2)*E_ij*DX(l,i,pi,nu)*DX(l,j,i,mu))/pow(B_ij,3) \
 						+ (4.0*E_ij*DX(l,i,pi,nu)*DX(l,j,i,mu))/pow(B_ij,2);
			
		}		
	}
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// ddVdr_nr
template <> void ddVdr_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = (j==0? (l.size()-1): j-1);
	uint pj = (j==(l.size()-1)? 0: j+1);		
	uint mk = (k==0? (l.size()-1): k-1);
	uint pk = (k==(l.size()-1)? 0: k+1);
	
	number B_jk = DistanceSquared(l[j],l[k]);
	number B_mjk = DistanceSquared(l[mj],l[k]);
	number B_jmk = DistanceSquared(l[j],l[mk]);
	number B_mjmk = DistanceSquared(l[mj],l[mk]);
	
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += + 2.0*pow(B_mjmk,-1.0 + a/2.0) \
					 + 2.0*pow(B_jk,-1.0 + a/2.0) \
					 + 4.0*pow(B_jk,-2.0 + a/2.0)*T_jk \
					 - 2.0*a*pow(B_jk,-2.0 + a/2.0)*T_jk ;
		if (k!=mj)
			res += - 2.0*pow(B_mjk,-1.0 + a/2.0); //
		if (k!=pj)
			res += - 2.0*pow(B_jmk,-1.0 + a/2.0); //
	}
	
	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res +=  - 4.0*pow(B_jk,-2.0 + a/2.0)*DX(l,j,pj,nu)*DX(l,j,k,mu) \
				 + 2.0*a*pow(B_jk,-2.0 + a/2.0)*DX(l,j,pj,nu)*DX(l,j,k,mu) \
				 + 4.0*pow(B_jk,-2.0 + a/2.0)*DX(l,j,k,nu)*DX(l,k,pk,mu) \
				 - 2.0*a*pow(B_jk,-2.0 + a/2.0)*DX(l,j,k,nu)*DX(l,k,pk,mu) \
				 - 16.0*pow(B_jk,-3.0 + a/2.0)*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk \
				 + 12.0*a*pow(B_jk,-3.0 + a/2.0)*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk \
				 - 2.0*pow(a,2)*pow(B_jk,-3.0 + a/2.0)*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk; //
	if (k!=mj)
		res += - 4.0*pow(B_mjk,-2.0 + a/2.0)*DX(l,mj,k,nu)*DX(l,k,pk,mu) \
 				+ 2.0*a*pow(B_mjk,-2.0 + a/2.0)*DX(l,mj,k,nu)*DX(l,k,pk,mu); //
	if (k!=pj)
		res += 	 + 4.0*pow(B_jmk,-2.0 + a/2.0)*DX(l,j,pj,nu)*DX(l,j,mk,mu) \
 				- 2.0*a*pow(B_jmk,-2.0 + a/2.0)*DX(l,j,pj,nu)*DX(l,j,mk,mu); //
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number B_ij, B_imj, T_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = (i==(l.size()-1)? 0: i+1);
			B_ij = DistanceSquared(l[i],l[j]);
			B_imj = DistanceSquared(l[i],l[mj]);
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);
			
			if (k==j && i!=j) {
				res += 	 - 4.0*pow(B_ij,-2.0 + a/2.0)*DX(l,i,pi,nu)*DX(l,j,i,mu) \
						 + 2.0*a*pow(B_ij,-2.0 + a/2.0)*DX(l,i,pi,nu)*DX(l,j,i,mu) \
						 - 4.0*pow(B_ij,-2.0 + a/2.0)*DX(l,i,pi,mu)*DX(l,j,i,nu) \
						 + 2.0*a*pow(B_ij,-2.0 + a/2.0)*DX(l,i,pi,mu)*DX(l,j,i,nu) \
						 + 16.0*pow(B_ij,-3.0 + a/2.0)*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij \
						 - 12.0*a*pow(B_ij,-3.0 + a/2.0)*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij \
						 + 2.0*pow(a,2)*pow(B_ij,-3.0 + a/2.0)*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij;
				if (mu==nu)
					res +=  - 4.0*pow(B_ij,-2.0 + a/2.0)*T_ij \
 							+ 2.0*a*pow(B_ij,-2.0 + a/2.0)*T_ij; //
			}
			if (k==mj && i!=mj) 
				res +=  4.0*pow(B_imj,-2.0 + a/2.0)*DX(l,i,pi,mu)*DX(l,mj,i,nu) \
 					- 2.0*a*pow(B_imj,-2.0 + a/2.0)*DX(l,i,pi,mu)*DX(l,mj,i,nu);
			if (k==pj && i!=j) 
				res += + 4.0*pow(B_ij,-2.0 + a/2.0)*DX(l,i,pi,nu)*DX(l,j,i,mu) \
 					- 2.0*a*pow(B_ij,-2.0 + a/2.0)*DX(l,i,pi,nu)*DX(l,j,i,mu);
			
		}		
	}
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// ddVthr_nr
template <> void ddVthr_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeigh(j,l.size());
	uint pj = posNeigh(j,l.size());		
	uint mk = negNeigh(k,l.size());
	uint pk = posNeigh(k,l.size());
	
	number r_jk = SpatialDistance(l[j],l[k]);
	number r_mjk = SpatialDistance(l[mj],l[k]);
	number r_jmk = SpatialDistance(l[j],l[mk]);
	number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
	number t_jk = DX(l,k,j,3);
	number t_mjk = DX(l,k,mj,3);
	number t_jmk = DX(l,mk,j,3);
	number t_mjmk = DX(l,mk,mj,3);
		
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*FThermal(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*FThermal(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*FThermal(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*FThermal(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DFThermalDrOnr_jk = DFThermalDrOnr(r_jk,t_jk,beta,a);
			number DFThermalDt_jk = DFThermalDt(r_jk,t_jk,beta,a);
			number DDFThermalDrDr_jk = DDFThermalDrDr(r_jk,t_jk,beta,a);
			number DDFThermalDrDtOnr_jk = DDFThermalDtDrOnr(r_jk,t_jk,beta,a);
			number DDFThermalDtDt_jk = DDFThermalDtDt(r_jk,t_jk,beta,a);
			if (mu==3) {
				res += - 2.0*DFThermalDt_jk*(-DX(l,j,nu)); //
				if (nu<3)
					res += (2.0*DDFThermalDrDtOnr_jk*DX(l,j,k,nu)*T_jk); //
				else
					res += -2.0*DDFThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFThermalDrOnr_jk*(-DX(l,j,nu))*DX(l,j,k,mu); //
				if (nu==3) 
					res += + 2.0*DDFThermalDrDtOnr_jk*DX(l,j,k,mu)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
					res += (2.0*DFThermalDrOnr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFThermalDrDr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DFThermalDt_jk*(-DX(l,k,mu)); //
			}
			else {
				res += - 2.0*DFThermalDrOnr_jk*DX(l,j,k,nu)*(-DX(l,k,mu)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DFThermalDt(r_mjk,t_mjk,beta,a)*(-DX(l,k,mu)); //
			else
				res += 2.0*DFThermalDrOnr(r_mjk,t_mjk,beta,a)*DX(l,mj,k,nu)*(-DX(l,k,mu)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DFThermalDt(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu));//
			else
				res += -(2.0*DFThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu))*DX(l,j,mk,mu)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DFThermalDrOnr_ij;	
		number DFThermalDt_ij;	
		number DDFThermalDrDr_ij;	
		number DDFThermalDrDtOnr_ij;
		number DDFThermalDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeigh(i,l.size());
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DX(l,j,i,3);
			t_imj = DX(l,mj,i,3);

			if (k==j && i!=j) {
				DFThermalDrOnr_ij = DFThermalDrOnr(r_ij,t_ij,beta,a);
				DFThermalDt_ij = DFThermalDt(r_ij,t_ij,beta,a);
				DDFThermalDrDr_ij = DDFThermalDrDr(r_ij,t_ij,beta,a);
				DDFThermalDrDtOnr_ij = DDFThermalDtDrOnr(r_ij,t_ij,beta,a);
				DDFThermalDtDt_ij = DDFThermalDtDt(r_ij,t_ij,beta,a);
				if (mu==3) {
					res += 2.0*DFThermalDt_ij*(-DX(l,i,nu)); //
					if (nu==3) 
						res += 2.0*DDFThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDFThermalDrDtOnr_ij*DX(l,j,i,nu)*T_ij; //
				}
				else {
					res += 2.0*DFThermalDrOnr_ij*(-DX(l,i,nu))*DX(l,j,i,mu); //
					if (nu==3)
						res += 2.0*DDFThermalDrDtOnr_ij*DX(l,j,i,mu)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
						res += - (2.0*DFThermalDrOnr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDFThermalDrDr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DFThermalDrOnr_ij*T_ij; //
				}
				if (nu==3) {
					res += 2.0*DFThermalDt_ij*(-DX(l,i,mu)); //
				}
				else {
					res += - 2.0*DFThermalDrOnr_ij*(-DX(l,i,mu))*DX(l,i,j,nu); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DFThermalDt(r_imj,t_imj,beta,a)*(-DX(l,i,mu));//
				else
					res += 2.0*DFThermalDrOnr(r_imj,t_imj,beta,a)*(-DX(l,i,mu))*DX(l,i,k,nu);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DFThermalDt(r_ij,t_ij,beta,a)*(-DX(l,i,nu));//
				else
					res += - 2.0*DFThermalDrOnr(r_ij,t_ij,beta,a)*(-DX(l,i,nu))*DX(l,j,i,mu);//
			}	
		}		
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (k==j && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*2.0/a/a;
	if (k==mj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	if (k==pj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	
	m(4*j+mu,4*k+nu) += f*(-pow(2.0*PI,2))*res;
}

// ddVthrDisjoint_nr
template <> void ddVthrDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeighDisjoint(j,l.size());
	uint pj = posNeighDisjoint(j,l.size());		
	uint mk = negNeighDisjoint(k,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	
	number r_jk = SpatialDistance(l[j],l[k]);
	number r_mjk = SpatialDistance(l[mj],l[k]);
	number r_jmk = SpatialDistance(l[j],l[mk]);
	number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
	number t_jk = DXDisjoint(l,k,j,3,beta);
	number t_mjk = DXDisjoint(l,k,mj,3,beta);
	number t_jmk = DXDisjoint(l,mk,j,3,beta);
	number t_mjmk = DXDisjoint(l,mk,mj,3,beta);
		
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*FThermal(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*FThermal(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*FThermal(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*FThermal(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DFThermalDrOnr_jk = DFThermalDrOnr(r_jk,t_jk,beta,a);
			number DFThermalDt_jk = DFThermalDt(r_jk,t_jk,beta,a);
			number DDFThermalDrDr_jk = DDFThermalDrDr(r_jk,t_jk,beta,a);
			number DDFThermalDrDtOnr_jk = DDFThermalDtDrOnr(r_jk,t_jk,beta,a);
			number DDFThermalDtDt_jk = DDFThermalDtDt(r_jk,t_jk,beta,a);
			if (mu==3) {
				res += - 2.0*DFThermalDt_jk*(-DXDisjoint(l,j,nu,beta)); //
				if (nu<3)
					res += (2.0*DDFThermalDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
				else
					res += -2.0*DDFThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFThermalDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
				if (nu==3) 
					res += + 2.0*DDFThermalDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
					res += (2.0*DFThermalDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFThermalDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DFThermalDt_jk*(-DXDisjoint(l,k,mu,beta)); //
			}
			else {
				res += - 2.0*DFThermalDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DFThermalDt(r_mjk,t_mjk,beta,a)*(-DXDisjoint(l,k,mu,beta)); //
			else
				res += 2.0*DFThermalDrOnr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DFThermalDt(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta));//
			else
				res += -(2.0*DFThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DFThermalDrOnr_ij;	
		number DFThermalDt_ij;	
		number DDFThermalDrDr_ij;	
		number DDFThermalDrDtOnr_ij;
		number DDFThermalDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeighDisjoint(i,l.size());
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DXDisjoint(l,j,i,3,beta);
			t_imj = DXDisjoint(l,mj,i,3,beta);

			if (k==j && i!=j) {
				DFThermalDrOnr_ij = DFThermalDrOnr(r_ij,t_ij,beta,a);
				DFThermalDt_ij = DFThermalDt(r_ij,t_ij,beta,a);
				DDFThermalDrDr_ij = DDFThermalDrDr(r_ij,t_ij,beta,a);
				DDFThermalDrDtOnr_ij = DDFThermalDtDrOnr(r_ij,t_ij,beta,a);
				DDFThermalDtDt_ij = DDFThermalDtDt(r_ij,t_ij,beta,a);
				if (mu==3) {
					res += 2.0*DFThermalDt_ij*(-DXDisjoint(l,i,nu,beta)); //
					if (nu==3) 
						res += 2.0*DDFThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDFThermalDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
				}
				else {
					res += 2.0*DFThermalDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
					if (nu==3)
						res += 2.0*DDFThermalDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
						res += - (2.0*DFThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDFThermalDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DFThermalDrOnr_ij*T_ij; //
				}
				if (nu==3) {
					res += 2.0*DFThermalDt_ij*(-DXDisjoint(l,i,mu,beta)); //
				}
				else {
					res += - 2.0*DFThermalDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DFThermalDt(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta));//
				else
					res += 2.0*DFThermalDrOnr(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DFThermalDt(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta));//
				else
					res += - 2.0*DFThermalDrOnr(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
			}	
		}		
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (k==j && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*2.0/a/a;
	if (k==mj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	if (k==pj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	
	m(4*j+mu,4*k+nu) += f*(-pow(2.0*PI,2))*res;
}

// ddVthrDisjointLR_nr
template <> void ddVthrDisjointLR_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	uint N = l.size();

	number a_jk = (((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))? a:0.0);
		
	number res = 0.0;
	
	uint mj = negNeighDisjoint(j,l.size());
	uint pj = posNeighDisjoint(j,l.size());		
	uint mk = negNeighDisjoint(k,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	
	number r_jk = SpatialDistance(l[j],l[k]);
	number r_mjk = SpatialDistance(l[mj],l[k]);
	number r_jmk = SpatialDistance(l[j],l[mk]);
	number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
	number t_jk = DXDisjoint(l,k,j,3,beta);
	number t_mjk = DXDisjoint(l,k,mj,3,beta);
	number t_jmk = DXDisjoint(l,mk,j,3,beta);
	number t_mjmk = DXDisjoint(l,mk,mj,3,beta);
		
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*FThermal(r_mjmk,t_mjmk,beta,a_jk)\
					+ 2.0*FThermal(r_jk,t_jk,beta,a_jk);//
		if (k!=mj)
			res +=  - 2.0*FThermal(r_mjk,t_mjk,beta,a_jk);//
		if (k!=pj)
			res += - 2.0*FThermal(r_jmk,t_jmk,beta,a_jk);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DFThermalDrOnr_jk = DFThermalDrOnr(r_jk,t_jk,beta,a_jk);
			number DFThermalDt_jk = DFThermalDt(r_jk,t_jk,beta,a_jk);
			number DDFThermalDrDr_jk = DDFThermalDrDr(r_jk,t_jk,beta,a_jk);
			number DDFThermalDrDtOnr_jk = DDFThermalDtDrOnr(r_jk,t_jk,beta,a_jk);
			number DDFThermalDtDt_jk = DDFThermalDtDt(r_jk,t_jk,beta,a_jk);
			if (mu==3) {
				res += - 2.0*DFThermalDt_jk*(-DXDisjoint(l,j,nu,beta)); //
				if (nu<3)
					res += (2.0*DDFThermalDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
				else
					res += -2.0*DDFThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFThermalDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
				if (nu==3) 
					res += + 2.0*DDFThermalDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
					res += (2.0*DFThermalDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFThermalDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DFThermalDt_jk*(-DXDisjoint(l,k,mu,beta)); //
			}
			else {
				res += - 2.0*DFThermalDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DFThermalDt(r_mjk,t_mjk,beta,a_jk)*(-DXDisjoint(l,k,mu,beta)); //
			else
				res += 2.0*DFThermalDrOnr(r_mjk,t_mjk,beta,a_jk)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DFThermalDt(r_jmk,t_jmk,beta,a_jk)*(-DXDisjoint(l,j,nu,beta));//
			else
				res += -(2.0*DFThermalDrOnr(r_jmk,t_jmk,beta,a_jk)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
			
		number a_ij;
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DFThermalDrOnr_ij;	
		number DFThermalDt_ij;	
		number DDFThermalDrDr_ij;	
		number DDFThermalDrDtOnr_ij;
		number DDFThermalDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeighDisjoint(i,l.size());
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DXDisjoint(l,j,i,3,beta);
			t_imj = DXDisjoint(l,mj,i,3,beta);
			a_ij = (((j<N/2 && i<N/2) || (j>=N/2 && i>=N/2))? a:0.0);

			if (k==j && i!=j) {
				DFThermalDrOnr_ij = DFThermalDrOnr(r_ij,t_ij,beta,a_ij);
				DFThermalDt_ij = DFThermalDt(r_ij,t_ij,beta,a_ij);
				DDFThermalDrDr_ij = DDFThermalDrDr(r_ij,t_ij,beta,a_ij);
				DDFThermalDrDtOnr_ij = DDFThermalDtDrOnr(r_ij,t_ij,beta,a_ij);
				DDFThermalDtDt_ij = DDFThermalDtDt(r_ij,t_ij,beta,a_ij);
				if (mu==3) {
					res += 2.0*DFThermalDt_ij*(-DXDisjoint(l,i,nu,beta)); //
					if (nu==3) 
						res += 2.0*DDFThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDFThermalDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
				}
				else {
					res += 2.0*DFThermalDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
					if (nu==3)
						res += 2.0*DDFThermalDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
						res += - (2.0*DFThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDFThermalDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DFThermalDrOnr_ij*T_ij; //
				}
				if (nu==3) {
					res += 2.0*DFThermalDt_ij*(-DXDisjoint(l,i,mu,beta)); //
				}
				else {
					res += - 2.0*DFThermalDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DFThermalDt(r_imj,t_imj,beta,a_ij)*(-DXDisjoint(l,i,mu,beta));//
				else
					res += 2.0*DFThermalDrOnr(r_imj,t_imj,beta,a_ij)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DFThermalDt(r_ij,t_ij,beta,a_ij)*(-DXDisjoint(l,i,nu,beta));//
				else
					res += - 2.0*DFThermalDrOnr(r_ij,t_ij,beta,a_ij)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
			}	
		}		
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (k==j && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*2.0/a/a;
	if (k==mj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	if (k==pj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	
	m(4*j+mu,4*k+nu) += f*(-pow(2.0*PI,2))*res;
}


// ddVnonrelDisjoint_nr
template <>
void ddVnonrelDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& f, mat& m) {
	
	if (mu!=3 && nu!=3) {
		uint oj = oppNeigh(j,l.size());
		uint ok = oppNeigh(k,l.size());
		number dt = 2.0*beta/(number)l.size();
	
		if (j==k) {
			m(4*j+mu,4*k+nu) += -f*dt*(-3.0)*DX(l,j,oj,mu)*DX(l,j,oj,nu)\
								/pow(SpatialDistance(l[j],l[oj]),5);
			if (mu==nu)
				m(4*j+mu,4*k+nu) += -f*dt/pow(SpatialDistance(l[j],l[oj]),3);
		}
		else if (k==oj) {
			m(4*j+mu,4*k+nu) += -f*dt*(-3.0)*DX(l,j,oj,mu)*DX(l,k,ok,nu)\
								/pow(SpatialDistance(l[j],l[oj]),5);
			if (mu==nu)
				m(4*j+mu,4*k+nu) += f*dt/pow(SpatialDistance(l[j],l[oj]),3);
		}
	}
	
}

// ddVnonrelrDisjoint_nr
template <>
void ddVnonrelrDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {

	if (mu!=3 && nu!=3) {
		uint oj = oppNeigh(j,l.size());
		uint ok = oppNeigh(k,l.size());
		number dt = 2.0*beta/(number)l.size();
		number r2 = SpatialDistanceSquared(l[j],l[oj]);
		number a2 = pow(a,2);	
	
		if (j==k) {
			m(4*j+mu,4*k+nu) += f*dt*DX(l,j,oj,mu)*DX(l,j,oj,nu)\
								*( 3.0/pow(r2+a2,5.0/2.0) - 4.0*exp(-r2/a2)/pow(a,5) );
			if (mu==nu)
				m(4*j+mu,4*k+nu) += f*dt*(-1.0/pow(r2 + a2,3.0/2.0) + 2.0*exp(-r2/a2)/pow(a,3) );
		}
		else if (k==oj) {
			m(4*j+mu,4*k+nu) += f*dt*DX(l,j,oj,mu)*DX(l,k,ok,nu)\
								*( 3.0/pow(r2+a2,5.0/2.0) - 4.0*exp(-r2/a2)/pow(a,5) );
			if (mu==nu)
				m(4*j+mu,4*k+nu) += -f*dt*(-1.0/pow(r2 + a2,3.0/2.0) + 2.0*exp(-r2/a2)/pow(a,3) );
		}
	}
	
}

// ddGaussian_nr
template <> void ddGaussian_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeigh(j,l.size());
	uint pj = posNeigh(j,l.size());		
	uint mk = negNeigh(k,l.size());
	uint pk = posNeigh(k,l.size());
	
	number B_jk = DistanceSquared(l[j],l[k]);
	number B_mjk = DistanceSquared(l[mj],l[k]);
	number B_jmk = DistanceSquared(l[j],l[mk]);
	number B_mjmk = DistanceSquared(l[mj],l[mk]);
	
	number E_jk = exp(-B_jk/a/a);
	number E_mjk = exp(-B_mjk/a/a);
	number E_jmk = exp(-B_jmk/a/a);
	number E_mjmk = exp(-B_mjmk/a/a);
	
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);

	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 	+ 2.0*E_mjmk \
					 + 2.0*E_jk \
					 + (4.0*E_jk*T_jk)/pow(a,2);
		if (k!=mj)
			res +=  - 2.0*E_mjk; //
		if (k!=pj)
			res += - 2.0*E_jmk; //
	}

	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res +=  - (4.0*E_jk*(-DX(l,j,nu))*DX(l,j,k,mu))/pow(a,2) \
				 + (4.0*E_jk*DX(l,j,k,nu)*(-DX(l,k,mu)))/pow(a,2) \
				 - (8.0*E_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(a,4); //
	if (k!=mj)
		res +=  - (4.0*E_mjk*DX(l,mj,k,nu)*(-DX(l,k,mu)))/pow(a,2); //
	if (k!=pj)
		res += + (4.0*E_jmk*(-DX(l,j,nu))*DX(l,j,mk,mu))/pow(a,2); //
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number B_ij, B_imj, E_ij, E_imj, T_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeigh(i,l.size());
			B_ij = DistanceSquared(l[i],l[j]);
			B_imj = DistanceSquared(l[i],l[mj]);
			E_ij = exp(-B_ij/a/a);
			E_imj = exp(-B_imj/a/a);
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);
			
			if (k==j && i!=j) {
				res += - (4.0*E_ij*(-DX(l,i,nu))*DX(l,j,i,mu))/pow(a,2) \
						 - (4.0*E_ij*(-DX(l,i,mu))*DX(l,j,i,nu))/pow(a,2) \
						 + (8.0*E_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(a,4);
				if (mu==nu)
					res +=  - (4.0*E_ij*T_ij)/pow(a,2); //
			}
			if (k==mj && i!=mj) 
				res +=  (4.0*E_imj*(-DX(l,i,mu))*DX(l,mj,i,nu))/pow(a,2);
			if (k==pj && i!=j) 
				res += + (4.0*E_ij*(-DX(l,i,nu))*DX(l,j,i,mu))/pow(a,2);
			
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0;
	if (k==mj && mu==nu)
		res += -1.0;
	if (k==pj && mu==nu)
		res += -1.0;
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// ddGaussianDisjoint_nr
template <> void ddGaussianDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeighDisjoint(j,l.size());
	uint pj = posNeighDisjoint(j,l.size());		
	uint mk = negNeighDisjoint(k,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	
	number B_jk = DistanceSquaredDisjoint(l[j],l[k],beta);
	number B_mjk = DistanceSquaredDisjoint(l[mj],l[k],beta);
	number B_jmk = DistanceSquaredDisjoint(l[j],l[mk],beta);
	number B_mjmk = DistanceSquaredDisjoint(l[mj],l[mk],beta);
	
	number E_jk = exp(-B_jk/a/a);
	number E_mjk = exp(-B_mjk/a/a);
	number E_jmk = exp(-B_jmk/a/a);
	number E_mjmk = exp(-B_mjmk/a/a);
	
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);

	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 	+ 2.0*E_mjmk \
					 + 2.0*E_jk \
					 + (4.0*E_jk*T_jk)/pow(a,2);
		if (k!=mj)
			res +=  - 2.0*E_mjk; //
		if (k!=pj)
			res += - 2.0*E_jmk; //
	}

	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res +=  - (4.0*E_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta))/pow(a,2) \
				 + (4.0*E_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)))/pow(a,2) \
				 - (8.0*E_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(a,4); //
	if (k!=mj)
		res +=  - (4.0*E_mjk*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)))/pow(a,2); //
	if (k!=pj)
		res += + (4.0*E_jmk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta))/pow(a,2); //
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number B_ij, B_imj, E_ij, E_imj, T_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeighDisjoint(i,l.size());
			B_ij = DistanceSquaredDisjoint(l[i],l[j],beta);
			B_imj = DistanceSquaredDisjoint(l[i],l[mj],beta);
			E_ij = exp(-B_ij/a/a);
			E_imj = exp(-B_imj/a/a);
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
			
			if (k==j && i!=j) {
				res += - (4.0*E_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta))/pow(a,2) \
						 - (4.0*E_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,j,i,nu,beta))/pow(a,2) \
						 + (8.0*E_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(a,4);
				if (mu==nu)
					res +=  - (4.0*E_ij*T_ij)/pow(a,2); //
			}
			if (k==mj && i!=mj) 
				res +=  (4.0*E_imj*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,mj,i,nu,beta))/pow(a,2);
			if (k==pj && i!=j) 
				res += + (4.0*E_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta))/pow(a,2);
			
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0;
	if (k==mj && mu==nu)
		res += -1.0;
	if (k==pj && mu==nu)
		res += -1.0;
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// ddGaussianLRDisjoint_nr
template <> void ddGaussianLRDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2)) {
		number res = 0.0;
	
		uint mj = negNeighDisjoint(j,l.size());
		uint pj = posNeighDisjoint(j,l.size());		
		uint mk = negNeighDisjoint(k,l.size());
		uint pk = posNeighDisjoint(k,l.size());
	
		number B_jk = DistanceSquaredDisjoint(l[j],l[k],beta);
		number B_mjk = DistanceSquaredDisjoint(l[mj],l[k],beta);
		number B_jmk = DistanceSquaredDisjoint(l[j],l[mk],beta);
		number B_mjmk = DistanceSquaredDisjoint(l[mj],l[mk],beta);
	
		number E_jk = exp(-B_jk/a/a);
		number E_mjk = exp(-B_mjk/a/a);
		number E_jmk = exp(-B_jmk/a/a);
		number E_mjmk = exp(-B_mjmk/a/a);
	
		number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);

		// terms where mu==nu, without sums
		if (mu==nu) {
			if (k!=j)
				res += 	+ 2.0*E_mjmk \
						 + 2.0*E_jk \
						 + (4.0*E_jk*T_jk)/pow(a,2);
			if (k!=mj)
				res +=  - 2.0*E_mjk; //
			if (k!=pj)
				res += - 2.0*E_jmk; //
		}

		// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j)
			res +=  - (4.0*E_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta))/pow(a,2) \
					 + (4.0*E_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)))/pow(a,2) \
					 - (8.0*E_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(a,4); //
		if (k!=mj)
			res +=  - (4.0*E_mjk*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)))/pow(a,2); //
		if (k!=pj)
			res += + (4.0*E_jmk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta))/pow(a,2); //
	
		// terms with sums
		if (k==j || k==mj || k==pj) {
	
			uint pi;
			number B_ij, B_imj, E_ij, E_imj, T_ij;
			
			// for LR
			uint iMin = (j<N/2? 0: N/2);
			uint iMax = (j<N/2? N/2:N);
		
			for (uint i=iMin; i<iMax; i++) {
		
				pi = posNeighDisjoint(i,l.size());
				B_ij = DistanceSquaredDisjoint(l[i],l[j],beta);
				B_imj = DistanceSquaredDisjoint(l[i],l[mj],beta);
				E_ij = exp(-B_ij/a/a);
				E_imj = exp(-B_imj/a/a);
				T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
			
				if (k==j && i!=j) {
					res += - (4.0*E_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta))/pow(a,2) \
							 - (4.0*E_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,j,i,nu,beta))/pow(a,2) \
							 + (8.0*E_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(a,4);
					if (mu==nu)
						res +=  - (4.0*E_ij*T_ij)/pow(a,2); //
				}
				if (k==mj && i!=mj) 
					res +=  (4.0*E_imj*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,mj,i,nu,beta))/pow(a,2);
				if (k==pj && i!=j) 
					res += + (4.0*E_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta))/pow(a,2);
			
			}		
		}
	
		//coincident terms
		if (k==j && mu==nu)
			res += 2.0;
		if (k==mj && mu==nu)
			res += -1.0;
		if (k==pj && mu==nu)
			res += -1.0;
	
		m(4*j+mu,4*k+nu) += f*res;
	}
}

// ddGaussianThermal_nr
template <> void ddGaussianThermal_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeigh(j,l.size());
	uint pj = posNeigh(j,l.size());		
	uint mk = negNeigh(k,l.size());
	uint pk = posNeigh(k,l.size());
	
	number r_jk = SpatialDistance(l[j],l[k]);
	number r_mjk = SpatialDistance(l[mj],l[k]);
	number r_jmk = SpatialDistance(l[j],l[mk]);
	number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
	number t_jk = DX(l,k,j,3);
	number t_mjk = DX(l,k,mj,3);
	number t_jmk = DX(l,mk,j,3);
	number t_mjmk = DX(l,mk,mj,3);
		
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*GThermal(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*GThermal(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*GThermal(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*GThermal(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DGThermalDrOnr_jk = DGThermalDrOnr(r_jk,t_jk,beta,a);
			number DGThermalDt_jk = DGThermalDt(r_jk,t_jk,beta,a);
			number DDGThermalDrDr_jk = DDGThermalDrDr(r_jk,t_jk,beta,a);
			number DDGThermalDrDtOnr_jk = DDGThermalDtDrOnr(r_jk,t_jk,beta,a);
			number DDGThermalDtDt_jk = DDGThermalDtDt(r_jk,t_jk,beta,a);
			if (mu==3) {
				res += - 2.0*DGThermalDt_jk*(-DX(l,j,nu)); //
				if (nu<3)
					res += (2.0*DDGThermalDrDtOnr_jk*DX(l,j,k,nu)*T_jk); //
				else
					res += -2.0*DDGThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DGThermalDrOnr_jk*(-DX(l,j,nu))*DX(l,j,k,mu); //
				if (nu==3) 
					res += + 2.0*DDGThermalDrDtOnr_jk*DX(l,j,k,mu)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DGThermalDrOnr+DDGThermalDrDr->0 as r->0
					res += (2.0*DGThermalDrOnr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2) \
							- (2.0*DDGThermalDrDr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DGThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DGThermalDt_jk*(-DX(l,k,mu)); //
			}
			else {
				res += - 2.0*DGThermalDrOnr_jk*DX(l,j,k,nu)*(-DX(l,k,mu)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DGThermalDt(r_mjk,t_mjk,beta,a)*(-DX(l,k,mu)); //
			else
				res += 2.0*DGThermalDrOnr(r_mjk,t_mjk,beta,a)*DX(l,mj,k,nu)*(-DX(l,k,mu)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DGThermalDt(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu));//
			else
				res += -(2.0*DGThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu))*DX(l,j,mk,mu)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DGThermalDrOnr_ij;	
		number DGThermalDt_ij;	
		number DDGThermalDrDr_ij;	
		number DDGThermalDrDtOnr_ij;
		number DDGThermalDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeigh(i,l.size());
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DX(l,j,i,3);
			t_imj = DX(l,mj,i,3);

			if (k==j && i!=j) {
				DGThermalDrOnr_ij = DGThermalDrOnr(r_ij,t_ij,beta,a);
				DGThermalDt_ij = DGThermalDt(r_ij,t_ij,beta,a);
				DDGThermalDrDr_ij = DDGThermalDrDr(r_ij,t_ij,beta,a);
				DDGThermalDrDtOnr_ij = DDGThermalDtDrOnr(r_ij,t_ij,beta,a);
				DDGThermalDtDt_ij = DDGThermalDtDt(r_ij,t_ij,beta,a);
				if (mu==3) {
					res += 2.0*DGThermalDt_ij*(-DX(l,i,nu)); //
					if (nu==3) 
						res += 2.0*DDGThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDGThermalDrDtOnr_ij*DX(l,j,i,nu)*T_ij; //
				}
				else {
					res += 2.0*DGThermalDrOnr_ij*(-DX(l,i,nu))*DX(l,j,i,mu); //
					if (nu==3)
						res += 2.0*DDGThermalDrDtOnr_ij*DX(l,j,i,mu)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DGThermalDrOnr+DDGThermalDrDr->0 as r->0
						res += - (2.0*DGThermalDrOnr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDGThermalDrDr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DGThermalDrOnr_ij*T_ij; //
				}
				if (nu==3) {
					res += 2.0*DGThermalDt_ij*(-DX(l,i,mu)); //
				}
				else {
					res += - 2.0*DGThermalDrOnr_ij*(-DX(l,i,mu))*DX(l,i,j,nu); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DGThermalDt(r_imj,t_imj,beta,a)*(-DX(l,i,mu));//
				else
					res += 2.0*DGThermalDrOnr(r_imj,t_imj,beta,a)*(-DX(l,i,mu))*DX(l,i,k,nu);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DGThermalDt(r_ij,t_ij,beta,a)*(-DX(l,i,nu));//
				else
					res += - 2.0*DGThermalDrOnr(r_ij,t_ij,beta,a)*(-DX(l,i,nu))*DX(l,j,i,mu);//
			}	
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0;
	if (k==mj && mu==nu)
		res += -1.0;
	if (k==pj && mu==nu)
		res += -1.0;
	
	m(4*j+mu,4*k+nu) += f*res;
}


// ddGaussianThermal2_nr
template <> void ddGaussianThermal2_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeigh(j,l.size());
	uint pj = posNeigh(j,l.size());		
	uint mk = negNeigh(k,l.size());
	uint pk = posNeigh(k,l.size());
	
	number r_jk = SpatialDistance(l[j],l[k]);
	number r_mjk = SpatialDistance(l[mj],l[k]);
	number r_jmk = SpatialDistance(l[j],l[mk]);
	number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
	number t_jk = DX(l,k,j,3);
	number t_mjk = DX(l,k,mj,3);
	number t_jmk = DX(l,mk,j,3);
	number t_mjmk = DX(l,mk,mj,3);
		
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*HThermal(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*HThermal(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*HThermal(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*HThermal(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DHThermalDrOnr_jk = DHThermalDrOnr(r_jk,t_jk,beta,a);
			number DHThermalDt_jk = DHThermalDt(r_jk,t_jk,beta,a);
			number DDHThermalDrDr_jk = DDHThermalDrDr(r_jk,t_jk,beta,a);
			number DDHThermalDrDtOnr_jk = DDHThermalDtDrOnr(r_jk,t_jk,beta,a);
			number DDHThermalDtDt_jk = DDHThermalDtDt(r_jk,t_jk,beta,a);
			if (mu==3) {
				res += - 2.0*DHThermalDt_jk*(-DX(l,j,nu)); //
				if (nu<3)
					res += (2.0*DDHThermalDrDtOnr_jk*DX(l,j,k,nu)*T_jk); //
				else
					res += -2.0*DDHThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DHThermalDrOnr_jk*(-DX(l,j,nu))*DX(l,j,k,mu); //
				if (nu==3) 
					res += + 2.0*DDHThermalDrDtOnr_jk*DX(l,j,k,mu)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DHThermalDrOnr+DDHThermalDrDr->0 as r->0
					res += (2.0*DHThermalDrOnr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2) \
							- (2.0*DDHThermalDrDr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DHThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DHThermalDt_jk*(-DX(l,k,mu)); //
			}
			else {
				res += - 2.0*DHThermalDrOnr_jk*DX(l,j,k,nu)*(-DX(l,k,mu)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DHThermalDt(r_mjk,t_mjk,beta,a)*(-DX(l,k,mu)); //
			else
				res += 2.0*DHThermalDrOnr(r_mjk,t_mjk,beta,a)*DX(l,mj,k,nu)*(-DX(l,k,mu)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DHThermalDt(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu));//
			else
				res += -(2.0*DHThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu))*DX(l,j,mk,mu)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DHThermalDrOnr_ij;	
		number DHThermalDt_ij;	
		number DDHThermalDrDr_ij;	
		number DDHThermalDrDtOnr_ij;
		number DDHThermalDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeigh(i,l.size());
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DX(l,j,i,3);
			t_imj = DX(l,mj,i,3);

			if (k==j && i!=j) {
				DHThermalDrOnr_ij = DHThermalDrOnr(r_ij,t_ij,beta,a);
				DHThermalDt_ij = DHThermalDt(r_ij,t_ij,beta,a);
				DDHThermalDrDr_ij = DDHThermalDrDr(r_ij,t_ij,beta,a);
				DDHThermalDrDtOnr_ij = DDHThermalDtDrOnr(r_ij,t_ij,beta,a);
				DDHThermalDtDt_ij = DDHThermalDtDt(r_ij,t_ij,beta,a);
				if (mu==3) {
					res += 2.0*DHThermalDt_ij*(-DX(l,i,nu)); //
					if (nu==3) 
						res += 2.0*DDHThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDHThermalDrDtOnr_ij*DX(l,j,i,nu)*T_ij; //
				}
				else {
					res += 2.0*DHThermalDrOnr_ij*(-DX(l,i,nu))*DX(l,j,i,mu); //
					if (nu==3)
						res += 2.0*DDHThermalDrDtOnr_ij*DX(l,j,i,mu)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DHThermalDrOnr+DDHThermalDrDr->0 as r->0
						res += - (2.0*DHThermalDrOnr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDHThermalDrDr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DHThermalDrOnr_ij*T_ij; //
				}
				if (nu==3) {
					res += 2.0*DHThermalDt_ij*(-DX(l,i,mu)); //
				}
				else {
					res += - 2.0*DHThermalDrOnr_ij*(-DX(l,i,mu))*DX(l,i,j,nu); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DHThermalDt(r_imj,t_imj,beta,a)*(-DX(l,i,mu));//
				else
					res += 2.0*DHThermalDrOnr(r_imj,t_imj,beta,a)*(-DX(l,i,mu))*DX(l,i,k,nu);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DHThermalDt(r_ij,t_ij,beta,a)*(-DX(l,i,nu));//
				else
					res += - 2.0*DHThermalDrOnr(r_ij,t_ij,beta,a)*(-DX(l,i,nu))*DX(l,j,i,mu);//
			}	
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0;
	if (k==mj && mu==nu)
		res += -1.0;
	if (k==pj && mu==nu)
		res += -1.0;
	
	m(4*j+mu,4*k+nu) += f*res;
}

// ddGaussianThermalDisjoint_nr
template <> void ddGaussianThermalDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeighDisjoint(j,l.size());
	uint pj = posNeighDisjoint(j,l.size());		
	uint mk = negNeighDisjoint(k,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	
	number r_jk = SpatialDistance(l[j],l[k]);
	number r_mjk = SpatialDistance(l[mj],l[k]);
	number r_jmk = SpatialDistance(l[j],l[mk]);
	number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
	number t_jk = DXDisjoint(l,k,j,3,beta);
	number t_mjk = DXDisjoint(l,k,mj,3,beta);
	number t_jmk = DXDisjoint(l,mk,j,3,beta);
	number t_mjmk = DXDisjoint(l,mk,mj,3,beta);
		
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*GThermal(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*GThermal(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*GThermal(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*GThermal(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DGThermalDrOnr_jk = DGThermalDrOnr(r_jk,t_jk,beta,a);
			number DGThermalDt_jk = DGThermalDt(r_jk,t_jk,beta,a);
			number DDGThermalDrDr_jk = DDGThermalDrDr(r_jk,t_jk,beta,a);
			number DDGThermalDrDtOnr_jk = DDGThermalDtDrOnr(r_jk,t_jk,beta,a);
			number DDGThermalDtDt_jk = DDGThermalDtDt(r_jk,t_jk,beta,a);
			if (mu==3) {
				res += - 2.0*DGThermalDt_jk*(-DXDisjoint(l,j,nu,beta)); //
				if (nu<3)
					res += (2.0*DDGThermalDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
				else
					res += -2.0*DDGThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DGThermalDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
				if (nu==3) 
					res += + 2.0*DDGThermalDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DGThermalDrOnr+DDGThermalDrDr->0 as r->0
					res += (2.0*DGThermalDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
							- (2.0*DDGThermalDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DGThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DGThermalDt_jk*(-DXDisjoint(l,k,mu,beta)); //
			}
			else {
				res += - 2.0*DGThermalDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DGThermalDt(r_mjk,t_mjk,beta,a)*(-DXDisjoint(l,k,mu,beta)); //
			else
				res += 2.0*DGThermalDrOnr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DGThermalDt(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta));//
			else
				res += -(2.0*DGThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DGThermalDrOnr_ij;	
		number DGThermalDt_ij;	
		number DDGThermalDrDr_ij;	
		number DDGThermalDrDtOnr_ij;
		number DDGThermalDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeighDisjoint(i,l.size());
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DXDisjoint(l,j,i,3,beta);
			t_imj = DXDisjoint(l,mj,i,3,beta);

			if (k==j && i!=j) {
				DGThermalDrOnr_ij = DGThermalDrOnr(r_ij,t_ij,beta,a);
				DGThermalDt_ij = DGThermalDt(r_ij,t_ij,beta,a);
				DDGThermalDrDr_ij = DDGThermalDrDr(r_ij,t_ij,beta,a);
				DDGThermalDrDtOnr_ij = DDGThermalDtDrOnr(r_ij,t_ij,beta,a);
				DDGThermalDtDt_ij = DDGThermalDtDt(r_ij,t_ij,beta,a);
				if (mu==3) {
					res += 2.0*DGThermalDt_ij*(-DXDisjoint(l,i,nu,beta)); //
					if (nu==3) 
						res += 2.0*DDGThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDGThermalDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
				}
				else {
					res += 2.0*DGThermalDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
					if (nu==3)
						res += 2.0*DDGThermalDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DGThermalDrOnr+DDGThermalDrDr->0 as r->0
						res += - (2.0*DGThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDGThermalDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DGThermalDrOnr_ij*T_ij; //
				}
				if (nu==3) {
					res += 2.0*DGThermalDt_ij*(-DXDisjoint(l,i,mu,beta)); //
				}
				else {
					res += - 2.0*DGThermalDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DGThermalDt(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta));//
				else
					res += 2.0*DGThermalDrOnr(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DGThermalDt(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta));//
				else
					res += - 2.0*DGThermalDrOnr(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
			}	
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0;
	if (k==mj && mu==nu)
		res += -1.0;
	if (k==pj && mu==nu)
		res += -1.0;
	
	m(4*j+mu,4*k+nu) += f*res;
}


// ddGaussianThermal2Disjoint_nr
template <> void ddGaussianThermal2Disjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeighDisjoint(j,l.size());
	uint pj = posNeighDisjoint(j,l.size());		
	uint mk = negNeighDisjoint(k,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	
	number r_jk = SpatialDistance(l[j],l[k]);
	number r_mjk = SpatialDistance(l[mj],l[k]);
	number r_jmk = SpatialDistance(l[j],l[mk]);
	number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
	number t_jk = DXDisjoint(l,k,j,3,beta);
	number t_mjk = DXDisjoint(l,k,mj,3,beta);
	number t_jmk = DXDisjoint(l,mk,j,3,beta);
	number t_mjmk = DXDisjoint(l,mk,mj,3,beta);
		
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*HThermal(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*HThermal(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*HThermal(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*HThermal(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DHThermalDrOnr_jk = DHThermalDrOnr(r_jk,t_jk,beta,a);
			number DHThermalDt_jk = DHThermalDt(r_jk,t_jk,beta,a);
			number DDHThermalDrDr_jk = DDHThermalDrDr(r_jk,t_jk,beta,a);
			number DDHThermalDrDtOnr_jk = DDHThermalDtDrOnr(r_jk,t_jk,beta,a);
			number DDHThermalDtDt_jk = DDHThermalDtDt(r_jk,t_jk,beta,a);
			if (mu==3) {
				res += - 2.0*DHThermalDt_jk*(-DXDisjoint(l,j,nu,beta)); //
				if (nu<3)
					res += (2.0*DDHThermalDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
				else
					res += -2.0*DDHThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DHThermalDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
				if (nu==3) 
					res += + 2.0*DDHThermalDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DHThermalDrOnr+DDHThermalDrDr->0 as r->0
					res += (2.0*DHThermalDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
							- (2.0*DDHThermalDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DHThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DHThermalDt_jk*(-DXDisjoint(l,k,mu,beta)); //
			}
			else {
				res += - 2.0*DHThermalDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DHThermalDt(r_mjk,t_mjk,beta,a)*(-DXDisjoint(l,k,mu,beta)); //
			else
				res += 2.0*DHThermalDrOnr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DHThermalDt(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta));//
			else
				res += -(2.0*DHThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DHThermalDrOnr_ij;	
		number DHThermalDt_ij;	
		number DDHThermalDrDr_ij;	
		number DDHThermalDrDtOnr_ij;
		number DDHThermalDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeighDisjoint(i,l.size());
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DXDisjoint(l,j,i,3,beta);
			t_imj = DXDisjoint(l,mj,i,3,beta);

			if (k==j && i!=j) {
				DHThermalDrOnr_ij = DHThermalDrOnr(r_ij,t_ij,beta,a);
				DHThermalDt_ij = DHThermalDt(r_ij,t_ij,beta,a);
				DDHThermalDrDr_ij = DDHThermalDrDr(r_ij,t_ij,beta,a);
				DDHThermalDrDtOnr_ij = DDHThermalDtDrOnr(r_ij,t_ij,beta,a);
				DDHThermalDtDt_ij = DDHThermalDtDt(r_ij,t_ij,beta,a);
				if (mu==3) {
					res += 2.0*DHThermalDt_ij*(-DXDisjoint(l,i,nu,beta)); //
					if (nu==3) 
						res += 2.0*DDHThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDHThermalDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
				}
				else {
					res += 2.0*DHThermalDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
					if (nu==3)
						res += 2.0*DDHThermalDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DHThermalDrOnr+DDHThermalDrDr->0 as r->0
						res += - (2.0*DHThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDHThermalDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DHThermalDrOnr_ij*T_ij; //
				}
				if (nu==3) {
					res += 2.0*DHThermalDt_ij*(-DXDisjoint(l,i,mu,beta)); //
				}
				else {
					res += - 2.0*DHThermalDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DHThermalDt(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta));//
				else
					res += 2.0*DHThermalDrOnr(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DHThermalDt(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta));//
				else
					res += - 2.0*DHThermalDrOnr(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
			}	
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += 2.0;
	if (k==mj && mu==nu)
		res += -1.0;
	if (k==pj && mu==nu)
		res += -1.0;
	
	m(4*j+mu,4*k+nu) += f*res;
}

// ddGaussianThermalLRDisjoint_nr
template <> void ddGaussianThermalLRDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2)) {
		number res = 0.0;
	
		uint mj = negNeighDisjoint(j,l.size());
		uint pj = posNeighDisjoint(j,l.size());		
		uint mk = negNeighDisjoint(k,l.size());
		uint pk = posNeighDisjoint(k,l.size());
	
		number r_jk = SpatialDistance(l[j],l[k]);
		number r_mjk = SpatialDistance(l[mj],l[k]);
		number r_jmk = SpatialDistance(l[j],l[mk]);
		number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
		number t_jk = DXDisjoint(l,k,j,3,beta);
		number t_mjk = DXDisjoint(l,k,mj,3,beta);
		number t_jmk = DXDisjoint(l,mk,j,3,beta);
		number t_mjmk = DXDisjoint(l,mk,mj,3,beta);
		
		number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
		// terms where mu==nu, without sums
		if (mu==nu) {
			if (k!=j)
				res += 2.0*GThermal(r_mjmk,t_mjmk,beta,a)\
						+ 2.0*GThermal(r_jk,t_jk,beta,a);//
			if (k!=mj)
				res +=  - 2.0*GThermal(r_mjk,t_mjk,beta,a);//
			if (k!=pj)
				res += - 2.0*GThermal(r_jmk,t_jmk,beta,a);//
		}

		// terms where mu not nexcessarily equal to nu, without sums
			if (k!=j) {
				number DGThermalDrOnr_jk = DGThermalDrOnr(r_jk,t_jk,beta,a);
				number DGThermalDt_jk = DGThermalDt(r_jk,t_jk,beta,a);
				number DDGThermalDrDr_jk = DDGThermalDrDr(r_jk,t_jk,beta,a);
				number DDGThermalDrDtOnr_jk = DDGThermalDtDrOnr(r_jk,t_jk,beta,a);
				number DDGThermalDtDt_jk = DDGThermalDtDt(r_jk,t_jk,beta,a);
				if (mu==3) {
					res += - 2.0*DGThermalDt_jk*(-DXDisjoint(l,j,nu,beta)); //
					if (nu<3)
						res += (2.0*DDGThermalDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
					else
						res += -2.0*DDGThermalDtDt_jk*T_jk;//
				}
				else {
					res += + 2.0*DGThermalDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
					if (nu==3) 
						res += + 2.0*DDGThermalDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
					else if (r_jk>MIN_NUMBER) //as -DGThermalDrOnr+DDGThermalDrDr->0 as r->0
						res += (2.0*DGThermalDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
								- (2.0*DDGThermalDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
					if (nu==mu)
						res += - (2.0*DGThermalDrOnr_jk*T_jk);//
				}
			
				if (nu==3) {
					res += 2.0*DGThermalDt_jk*(-DXDisjoint(l,k,mu,beta)); //
				}
				else {
					res += - 2.0*DGThermalDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
				}
			
			}
			if (k!=mj) {
				if (nu==3)
					res += - 2.0*DGThermalDt(r_mjk,t_mjk,beta,a)*(-DXDisjoint(l,k,mu,beta)); //
				else
					res += 2.0*DGThermalDrOnr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			if (j!=mk) {
				if (mu==3)
					res += 2.0*DGThermalDt(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta));//
				else
					res += -(2.0*DGThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
			}
	
		// terms with sums
		if (k==j || k==mj || k==pj) {
	
			uint pi;
			number T_ij;
			number r_ij, r_imj;
			number t_ij, t_imj;
			number DGThermalDrOnr_ij;	
			number DGThermalDt_ij;	
			number DDGThermalDrDr_ij;	
			number DDGThermalDrDtOnr_ij;
			number DDGThermalDtDt_ij;
			
			// for LR
			uint iMin = (j<N/2? 0: N/2);
			uint iMax = (j<N/2? N/2: N);
		
			for (uint i=iMin; i<iMax; i++) {
		
				pi = posNeighDisjoint(i,l.size());
				T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
				r_ij = SpatialDistance(l[j],l[i]);
				r_imj = SpatialDistance(l[i],l[mj]);	
				t_ij = DXDisjoint(l,j,i,3,beta);
				t_imj = DXDisjoint(l,mj,i,3,beta);

				if (k==j && i!=j) {
					DGThermalDrOnr_ij = DGThermalDrOnr(r_ij,t_ij,beta,a);
					DGThermalDt_ij = DGThermalDt(r_ij,t_ij,beta,a);
					DDGThermalDrDr_ij = DDGThermalDrDr(r_ij,t_ij,beta,a);
					DDGThermalDrDtOnr_ij = DDGThermalDtDrOnr(r_ij,t_ij,beta,a);
					DDGThermalDtDt_ij = DDGThermalDtDt(r_ij,t_ij,beta,a);
					if (mu==3) {
						res += 2.0*DGThermalDt_ij*(-DXDisjoint(l,i,nu,beta)); //
						if (nu==3) 
							res += 2.0*DDGThermalDtDt_ij*T_ij; //
						else
							res += 2.0*DDGThermalDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
					}
					else {
						res += 2.0*DGThermalDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
						if (nu==3)
							res += 2.0*DDGThermalDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
						else if (r_ij>MIN_NUMBER) { // as -DGThermalDrOnr+DDGThermalDrDr->0 as r->0
							res += - (2.0*DGThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
								+ (2.0*DDGThermalDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
						}
						if (nu==mu)
							res += 2.0*DGThermalDrOnr_ij*T_ij; //
					}
					if (nu==3) {
						res += 2.0*DGThermalDt_ij*(-DXDisjoint(l,i,mu,beta)); //
					}
					else {
						res += - 2.0*DGThermalDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
					}
				} 
				if (j==pk && i!=mj) {
					if (nu==3)
						res += - 2.0*DGThermalDt(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta));//
					else
						res += 2.0*DGThermalDrOnr(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
				} 
				if (j==mk && i!=j) {
					if (mu==3)
						res += - 2.0*DGThermalDt(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta));//
					else
						res += - 2.0*DGThermalDrOnr(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
				}	
			}		
		}
	
		//coincident terms
		if (k==j && mu==nu)
			res += 2.0;
		if (k==mj && mu==nu)
			res += -1.0;
		if (k==pj && mu==nu)
			res += -1.0;
	
		m(4*j+mu,4*k+nu) += f*res;
	}
}


// ddGaussianThermal2LRDisjoint_nr
template <> void ddGaussianThermal2LRDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2)) {
		number res = 0.0;
	
		uint mj = negNeighDisjoint(j,l.size());
		uint pj = posNeighDisjoint(j,l.size());		
		uint mk = negNeighDisjoint(k,l.size());
		uint pk = posNeighDisjoint(k,l.size());
	
		number r_jk = SpatialDistance(l[j],l[k]);
		number r_mjk = SpatialDistance(l[mj],l[k]);
		number r_jmk = SpatialDistance(l[j],l[mk]);
		number r_mjmk = SpatialDistance(l[mj],l[mk]);
	
		number t_jk = DXDisjoint(l,k,j,3,beta);
		number t_mjk = DXDisjoint(l,k,mj,3,beta);
		number t_jmk = DXDisjoint(l,mk,j,3,beta);
		number t_mjmk = DXDisjoint(l,mk,mj,3,beta);
		
		number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
		// terms where mu==nu, without sums
		if (mu==nu) {
			if (k!=j)
				res += 2.0*HThermal(r_mjmk,t_mjmk,beta,a)\
						+ 2.0*HThermal(r_jk,t_jk,beta,a);//
			if (k!=mj)
				res +=  - 2.0*HThermal(r_mjk,t_mjk,beta,a);//
			if (k!=pj)
				res += - 2.0*HThermal(r_jmk,t_jmk,beta,a);//
		}

		// terms where mu not nexcessarily equal to nu, without sums
			if (k!=j) {
				number DHThermalDrOnr_jk = DHThermalDrOnr(r_jk,t_jk,beta,a);
				number DHThermalDt_jk = DHThermalDt(r_jk,t_jk,beta,a);
				number DDHThermalDrDr_jk = DDHThermalDrDr(r_jk,t_jk,beta,a);
				number DDHThermalDrDtOnr_jk = DDHThermalDtDrOnr(r_jk,t_jk,beta,a);
				number DDHThermalDtDt_jk = DDHThermalDtDt(r_jk,t_jk,beta,a);
				if (mu==3) {
					res += - 2.0*DHThermalDt_jk*(-DXDisjoint(l,j,nu,beta)); //
					if (nu<3)
						res += (2.0*DDHThermalDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
					else
						res += -2.0*DDHThermalDtDt_jk*T_jk;//
				}
				else {
					res += + 2.0*DHThermalDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
					if (nu==3) 
						res += + 2.0*DDHThermalDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
					else if (r_jk>MIN_NUMBER) //as -DHThermalDrOnr+DDHThermalDrDr->0 as r->0
						res += (2.0*DHThermalDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
								- (2.0*DDHThermalDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
					if (nu==mu)
						res += - (2.0*DHThermalDrOnr_jk*T_jk);//
				}
			
				if (nu==3) {
					res += 2.0*DHThermalDt_jk*(-DXDisjoint(l,k,mu,beta)); //
				}
				else {
					res += - 2.0*DHThermalDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
				}
			
			}
			if (k!=mj) {
				if (nu==3)
					res += - 2.0*DHThermalDt(r_mjk,t_mjk,beta,a)*(-DXDisjoint(l,k,mu,beta)); //
				else
					res += 2.0*DHThermalDrOnr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			if (j!=mk) {
				if (mu==3)
					res += 2.0*DHThermalDt(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta));//
				else
					res += -(2.0*DHThermalDrOnr(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
			}
	
		// terms with sums
		if (k==j || k==mj || k==pj) {
	
			uint pi;
			number T_ij;
			number r_ij, r_imj;
			number t_ij, t_imj;
			number DHThermalDrOnr_ij;	
			number DHThermalDt_ij;	
			number DDHThermalDrDr_ij;	
			number DDHThermalDrDtOnr_ij;
			number DDHThermalDtDt_ij;
			
			// for LR
			uint iMin = (j<N/2? 0: N/2);
			uint iMax = (j<N/2? N/2: N);
		
			for (uint i=iMin; i<iMax; i++) {
		
				pi = posNeighDisjoint(i,l.size());
				T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
				r_ij = SpatialDistance(l[j],l[i]);
				r_imj = SpatialDistance(l[i],l[mj]);	
				t_ij = DXDisjoint(l,j,i,3,beta);
				t_imj = DXDisjoint(l,mj,i,3,beta);

				if (k==j && i!=j) {
					DHThermalDrOnr_ij = DHThermalDrOnr(r_ij,t_ij,beta,a);
					DHThermalDt_ij = DHThermalDt(r_ij,t_ij,beta,a);
					DDHThermalDrDr_ij = DDHThermalDrDr(r_ij,t_ij,beta,a);
					DDHThermalDrDtOnr_ij = DDHThermalDtDrOnr(r_ij,t_ij,beta,a);
					DDHThermalDtDt_ij = DDHThermalDtDt(r_ij,t_ij,beta,a);
					if (mu==3) {
						res += 2.0*DHThermalDt_ij*(-DXDisjoint(l,i,nu,beta)); //
						if (nu==3) 
							res += 2.0*DDHThermalDtDt_ij*T_ij; //
						else
							res += 2.0*DDHThermalDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
					}
					else {
						res += 2.0*DHThermalDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
						if (nu==3)
							res += 2.0*DDHThermalDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
						else if (r_ij>MIN_NUMBER) { // as -DHThermalDrOnr+DDHThermalDrDr->0 as r->0
							res += - (2.0*DHThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
								+ (2.0*DDHThermalDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
						}
						if (nu==mu)
							res += 2.0*DHThermalDrOnr_ij*T_ij; //
					}
					if (nu==3) {
						res += 2.0*DHThermalDt_ij*(-DXDisjoint(l,i,mu,beta)); //
					}
					else {
						res += - 2.0*DHThermalDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
					}
				} 
				if (j==pk && i!=mj) {
					if (nu==3)
						res += - 2.0*DHThermalDt(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta));//
					else
						res += 2.0*DHThermalDrOnr(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
				} 
				if (j==mk && i!=j) {
					if (mu==3)
						res += - 2.0*DHThermalDt(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta));//
					else
						res += - 2.0*DHThermalDrOnr(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
				}	
			}		
		}
	
		//coincident terms
		if (k==j && mu==nu)
			res += 2.0;
		if (k==mj && mu==nu)
			res += -1.0;
		if (k==pj && mu==nu)
			res += -1.0;
	
		m(4*j+mu,4*k+nu) += f*res;
	}
}

