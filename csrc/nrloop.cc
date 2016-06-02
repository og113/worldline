/*
	definitions for solving classical equations of motion for monopole worldline via newton-raphson method
*/

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
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
		4 - explicit instatiation
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0 - some static functions
----------------------------------------------------------------------------------------------------------------------------*/

// sinch
static number sinch(const number& x) {
	return (abs(x)<MIN_NUMBER? 1.0: sinh(x)/x);
}

// coshsinhc
static number coshsinhc(const number& x) {
	return (abs(x)<MIN_NUMBER? 0.0: cosh(x)/x - sinh(x)/pow(x,2));
}

// coshsinhcc
static number coshsinhcc(const number& x) {
	return (abs(x)<MIN_NUMBER? 1.0/3.0: cosh(x)/pow(x,2) - sinh(x)/pow(x,3));
}

// FThermal
number FThermal(const number& r, const number& t, const number& beta, const number& a) {
	return sinch((2.0*r*PI)/beta)/(2.0*pow(beta,2)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*r*PI)/beta)));
}

// DFThermalDr
number DFThermalDr(const number& r, const number& t, const number& beta, const number& a) {
	return (PI*coshsinhc((2.0*PI*r)/beta))/\
 (pow(beta,3)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta))) \
 + (2.0*pow(PI,2)*r*pow(sinch((2.0*PI*r)/beta),2))/\
 (pow(beta,4)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
}

// DFThermalDrOnr
number DFThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (2.0*pow(PI,2)*coshsinhcc((2.0*PI*r)/beta))/\
 (pow(beta,4)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta))) \
 + (2.0*pow(PI,2)*pow(sinch((2.0*PI*r)/beta),2))/\
 (pow(beta,4)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
 }

// DFThermalDt
number DFThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	return (PI*sin((2.0*PI*t)/beta)*sinch((2.0*PI*r)/beta))/\
 (pow(beta,3)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
}

// DDFThermalDrDr
number DDFThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*pow(PI,2)*(-((2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta))*
 (3.0*cosh((2.0*PI*r)/beta) - 2.0*sinch((2.0*PI*r)/beta))*sinch((2.0*PI*r)/beta)) \
 + 8.0*pow(PI,2)*pow(r,2)*pow(sinch((2.0*PI*r)/beta),3) \
 + (pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),2)*
 (-2.0*coshsinhcc((2.0*PI*r)/beta) + sinch((2.0*PI*r)/beta)))/pow(beta,2)))/\
 pow(-2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cos((2.0*PI*t)/beta) - pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}

// DDFThermalDtDr
number DDFThermalDtDr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*pow(PI,2)*sin((2.0*PI*t)/beta)*((2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta))*
 coshsinhc((2.0*PI*r)/beta) - 4.0*PI*r*beta*pow(sinch((2.0*PI*r)/beta),2)))/\
 pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}

 // DDFThermalDtDrOnr
number DDFThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (4.0*pow(PI,3)*sin((2.0*PI*t)/beta)*((-2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cos((2.0*PI*t)/beta) - pow(beta,2)*cosh((2.0*PI*r)/beta))*\
 coshsinhcc((2.0*PI*r)/beta) + 2.0*pow(beta,2)*pow(sinch((2.0*PI*r)/beta),2)))/\
 (beta*pow(-2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cos((2.0*PI*t)/beta) - pow(beta,2)*cosh((2.0*PI*r)/beta),3));
}

// DDFThermalDtDt
number DDFThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	return  (pow(PI,2)*(-3.0*pow(beta,2) + pow(beta,2)*cos((4.0*PI*t)/beta) \
 + 2.0*cos((2.0*PI*t)/beta)*(2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cosh((2.0*PI*r)/beta)))*sinch((2.0*PI*r)/beta))/\
 pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}

/*----------------------------------------------------------------------------------------------------------------------------
	1 - nr loop functions
----------------------------------------------------------------------------------------------------------------------------*/

// L
template <uint Dim>
void L (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pos = (j==(l.size()-1)? 0: j+1);
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
	uint pos = (j==(l.size()-1)? 0: j+1);
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
	uint pj = (j==(l.size()-1)? 0: j+1);
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
	cerr << "I0 error: not defined for dim = " << Dim << endl;
}

// I0Disjoint
template <uint Dim>
void I0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	cerr << "I0Disjoint error: not defined for dim = " << Dim << endl;
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
		cerr << "MaxXn error: " << n << ">" << Dim << endl;
	}
	if (k<j) {
		number temp = abs((l[j])[n]-(l[k])[n]);
		result = (temp*f>result? temp*f: result);
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
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
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
	
		result += f*(1.0+(number)(k<j))*DotDisjoint(l[pj],l[j],l[pk],l[k],beta)\
						*exp(-DistanceSquaredDisjoint(l[j],l[k],beta)/a/a);	
	}
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

// mdI_nr
template<uint Dim>
void mdI_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
	cerr << "mdI_nr error: not defined for dim = " << Dim << endl;
}

// ddI_nr
template<uint Dim>
void ddI_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
	cerr << "ddI_nr error: not defined for dim = " << Dim << endl;
}

// mdIDisjoint_nr
template<uint Dim>
void mdIDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	cerr << "mdIDisjoint_nr error: not defined for dim = " << Dim << endl;
}

// ddIDisjoint_nr
template<uint Dim>
void ddIDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& beta, const number& f, mat& m) {
	cerr << "ddIDisjoint_nr error: not defined for dim = " << Dim << endl;
}

// mdS0_nr
template<uint Dim>
void mdS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	v[j*Dim+mu] += -f*(DX(l,j,nj,mu) + DX(l,j,pj,mu))*(number)l.size()/2.0;
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

// mdsqrtS0_nr
template<uint Dim>
void mdsqrtS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0, const number& f, vec& v) {
	mdS0_nr(j,mu,l,f*2.0/sqrt4s0,v);
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

// mdS0Disjoint_nr
template<uint Dim>
void mdS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	v[j*Dim+mu] += -f*(DXDisjoint(l,j,nj,mu,beta) + DXDisjoint(l,j,pj,mu,beta))*(number)l.size()/2.0;
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

// mdsqrtS0Disjoint_nr
template<uint Dim>
void mdsqrtS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0,\
 const number& beta, const number& f, vec& v) {
 	mdS0Disjoint_nr(j,mu,l,beta,f*2.0/sqrt4s0,v);
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
	
	m(j*Dim+mu,k*Dim+nu) -= (f*pow((number)l.size(),2)/pow(sqrt4s0,3)) * \
									(DXDisjoint(l,j,nj,mu,beta) + DXDisjoint(l,j,pj,mu,beta)) \
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

// ddVthr_nr
template<uint Dim>
void ddVthr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& theta, const number& a, const number& f, mat& m) {
	cerr << "mddVthr_nr Error: no script written for dim = " << Dim << endl;
}

// mdVthrDisjoint_nr
template<uint Dim>
void mdVthr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mddVthrDisjoint_nr Error: no script written for dim = " << Dim << endl;
}

// ddVthrDisjoint_nr
template<uint Dim>
void ddVthrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "mddVthrDisjoint_nr Error: no script written for dim = " << Dim << endl;
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

// ddGaussian_nr
template<uint Dim>
void ddGaussian_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "ddGaussian_nr Error: no script written for dim = " << Dim << endl;
}

// mdGaussianDisjoint_nr
template<uint Dim>
void mdGaussianDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
		for (uint i=0; i<l.size(); i++)
	mdGaussian_nrDisjoint( j, mu, i, l, a, f, v);
}

// mdGaussianDisjoint_nr
template<uint Dim>
void mdGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	cerr << "mdGaussianDisjoint_nr Error: no script written for dim = " << Dim << endl;
}
						 
// ddGaussianDisjoint_nr
template<uint Dim>
void ddGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	cerr << "ddGaussianDisjoint_nr Error: no script written for dim = " << Dim << endl;
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

// PS0_nr
template<uint Dim>
void PS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v) {
	v[Dim*loc+mu] += f*DX(l,posNeigh(loc,l.size()),loc,mu)*(number)l.size()/2.0;
}

// PS0Disjoint_nr
template<uint Dim>
void PS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v) {
	v[Dim*loc+mu] += f*DXDisjoint(l,posNeighDisjoint(loc,l.size()),loc,mu,beta)*(number)l.size()/2.0;
}

// PL_nr
template<uint Dim>
void PL_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v) {
	uint ploc = posNeigh(loc,l.size());
	v[Dim*loc+mu] += f*DX(l,ploc,loc,mu)/Distance(l[ploc],l[loc]);
}

// PLDisjoint_nr
template<uint Dim>
void PLDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v) {
	uint ploc = posNeighDisjoint(loc,l.size());
	v[Dim*loc+mu] += f*DXDisjoint(l,ploc,loc,mu,beta)/DistanceDisjoint(l[ploc],l[loc],beta);
}

// PsqrtS0_nr
template<uint Dim>
void PsqrtS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, vec& v) {
	PS0_nr(l,loc,mu,f*2.0/sqrt4s0,v);
}

// PsqrtS0Disjoint_nr
template<uint Dim>
void PsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, vec& v) {
	PS0Disjoint_nr(l,loc,mu,beta,f*2.0/sqrt4s0,v);
}

// PI0_nr
template<uint Dim>
void PI0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v) {
	if (mu==(Dim-2)) 
		v[Dim*loc+mu] += - f*(l[loc])[Dim-1];
	else if (mu==(Dim-1))
		v[Dim*loc+mu] += f*(l[loc])[Dim-2];
}

// PI0Disjoint_nr
template<uint Dim>
void PI0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v) {
	PI0_nr(l,loc,mu,f,v);
}

// PVor_nr
template<uint Dim>
void PVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v) {
	uint pk = posNeigh(k,l.size());
	v[Dim*j+mu] += f*2.0*DX(l,pk,k)/(DistanceSquared(l[j],l[k]) + a*a);
}

// PVthr_nr
template<uint Dim>
void PVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint pj = posNeigh(j,l.size());
	uint pk = posNeigh(k,l.size());
	
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);

	v[Dim*j+mu] += f*(-pow(2.0*PI,2))*2.0*DX(l,pk,k)*FThermal(r,t,beta,a);	
}

// PVthrDisjoint_nr
template<uint Dim>
void PVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint pj = posNeighDisjoint(j,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	
	number r = SpatialDistance(l[j],l[k]);
	number t = DXDisjoint(l,k,j,Dim-1,beta);

	v[Dim*j+mu] += f*(-pow(2.0*PI,2))*2.0*DXDisjoint(l,pk,k,beta)*FThermal(r,t,beta,a);	
}

// PGaussian_nr
template<uint Dim>
void PGaussian_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v) {
	uint pk = posNeigh(k,l.size());
	v[Dim*j+mu] += f*2.0*DX(l,pk,k)*exp(-DistanceSquared(l[j],l[k])/a/a);
}

// PGaussianDisjoint_nr
template<uint Dim>
void PGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size());
	v[Dim*j+mu] += f*2.0*DXDisjoint(l,pk,k,beta)*exp(-DistanceSquaredDisjoint(l[j],l[k],beta)/a/a);
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
Filename filenameLoopNR(const Parameters& p) {
	Filename f = "data/nr/loops/dim_"+nts(Dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(p.P4)\
		+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+".dat";
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
template Filename filenameLoopNR<2>(const Parameters& p);

// mdI_nr<2>
template <> void mdI_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& f, vec& v) {
	if (mu==1) {
		uint nj = (j==0? (l.size()-1): j-1);
		v[j*2+mu] += -f*((l[nj])[0]-(l[j])[0]);
	}
	else if (mu==0) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		v[j*2+mu] += -f*((l[pj])[1]-(l[j])[1]);
	}
}

// ddI_nr<2>
template <> void ddI_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
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
template void Sm<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void Angle<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void FGamma<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void Gaussian<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result);
template void GaussianDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void MaxXn<4>(const uint& j, const uint& k, const Loop<4>& l, const uint& n, const number& f, number& result);
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
template void mdVor_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& a, const number& p, vec& v);
template void mdFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, vec& v);
template void ddFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, mat& v);
template void loopToVector<4>(const Loop<4>&,vec&);
template void vectorToLoop<4>(const vec&, Loop<4>&);
template Filename filenameLoopNR<4>(const Parameters& p);

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
		
		res += 2.0*FThermal_ij*DX(l,i,pi,mu);
		if (mu<3)
			res += 2.0*DFThermalDrOnr_ij*DX(l,j,i,mu)*T_ij;
		else
			res += 2.0*DFThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DX(l,mj,i,3); // checked order
		number FThermal_imj = FThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*FThermal_imj*DX(l,i,pi,mu); //
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
		
		res += 2.0*FThermal_ij*DXDisjoint(l,i,pi,mu,beta);
		if (mu<3)
			res += 2.0*DFThermalDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij;
		else
			res += 2.0*DFThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DXDisjoint(l,mj,i,3,beta); // checked order
		number FThermal_imj = FThermal(r_imj,t_imj,beta,a);
		res +=  -2.0*FThermal_imj*DXDisjoint(l,i,pi,mu,beta); //
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
		res += + 2.0*E_ij*DX(l,i,pi,mu) \
 			- (4.0*E_ij*DX(l,j,i,mu)*T_ij)/pow(a,2);
	}
	if (i!=mj) {
		number B_imj = DistanceSquared(l[i],l[mj]);
		number E_imj = exp(-B_imj/a/a);
		res += - 2.0*E_imj*DX(l,i,pi,mu);
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
		res += + 2.0*E_ij*DXDisjoint(l,i,pi,mu,beta) \
 			- (4.0*E_ij*DXDisjoint(l,j,i,mu,beta)*T_ij)/pow(a,2);
	}
	if (i!=mj) {
		number B_imj = DistanceSquaredDisjoint(l[i],l[mj],beta);
		number E_imj = exp(-B_imj/a/a);
		res += - 2.0*E_imj*DXDisjoint(l,i,pi,mu,beta);
	}
		
	//coincident terms	
	if (mu==3) {
		if (i==j) { // bit of a fudge in terms of sum over i
			res += (DXDisjoint(l,j,mj,mu,beta) + DXDisjoint(l,j,pj,mu,beta))/a/a;
		}
	}
	else {
		if (i==j)
			res += 2.0*(l[j])[mu]/a/a;
		if (i==mj)
			res += -(l[mj])[mu]/a/a;
		if (i==pj)
			res += -(l[pj])[mu]/a/a;
	}

	v[j*4+mu] += -f*res;	

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
				res += - 2.0*DFThermalDt_jk*DX(l,j,pj,nu); //
				if (nu<3)
					res += (2.0*DDFThermalDrDtOnr_jk*DX(l,j,k,nu)*T_jk); //
				else
					res += -2.0*DDFThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFThermalDrOnr_jk*DX(l,j,pj,nu)*DX(l,j,k,mu); //
				if (nu==3) 
					res += + 2.0*DDFThermalDrDtOnr_jk*DX(l,j,k,mu)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
					res += (2.0*DFThermalDrOnr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFThermalDrDr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DFThermalDt_jk*DX(l,k,pk,mu); //
			}
			else {
				res += - (2.0*DFThermalDrOnr_jk*DX(l,j,k,nu)*DX(l,k,pk,mu)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DFThermalDt(r_mjk,t_mjk,beta,a)*DX(l,k,pk,mu); //
			else
				res += (2.0*DFThermalDrOnr(r_mjk,t_mjk,beta,a)*DX(l,mj,k,nu)*DX(l,k,pk,mu)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DFThermalDt(r_jmk,t_jmk,beta,a)*DX(l,j,pj,nu);//
			else
				res += -(2.0*DFThermalDrOnr(r_jmk,t_jmk,beta,a)*DX(l,j,pj,nu)*DX(l,j,mk,mu)); //
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
					res += 2.0*DFThermalDt_ij*DX(l,i,pi,nu); //
					if (nu==3) 
						res += 2.0*DDFThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDFThermalDrDtOnr_ij*DX(l,j,i,nu)*T_ij; //
				}
				else {
					res += 2.0*DFThermalDrOnr_ij*DX(l,i,pi,nu)*DX(l,j,i,mu); //
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
					res += 2.0*DFThermalDt_ij*DX(l,i,pi,mu); //
				}
				else {
					res += - 2.0*DFThermalDrOnr_ij*DX(l,i,pi,mu)*DX(l,i,j,nu); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DFThermalDt(r_imj,t_imj,beta,a)*DX(l,i,pi,mu);//
				else
					res += 2.0*DFThermalDrOnr(r_imj,t_imj,beta,a)*DX(l,i,pi,mu)*DX(l,i,k,nu);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DFThermalDt(r_ij,t_ij,beta,a)*DX(l,i,pi,nu);//
				else
					res += - 2.0*DFThermalDrOnr(r_ij,t_ij,beta,a)*DX(l,i,pi,nu)*DX(l,j,i,mu);//
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
				res += - 2.0*DFThermalDt_jk*DXDisjoint(l,j,pj,nu,beta); //
				if (nu<3)
					res += (2.0*DDFThermalDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
				else
					res += -2.0*DDFThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFThermalDrOnr_jk*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,k,mu,beta); //
				if (nu==3) 
					res += + 2.0*DDFThermalDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
					res += (2.0*DFThermalDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFThermalDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==3) {
				res += 2.0*DFThermalDt_jk*DXDisjoint(l,k,pk,mu,beta); //
			}
			else {
				res += - (2.0*DFThermalDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta)); //
			}
			
		}
		if (k!=mj) {
			if (nu==3)
				res += - 2.0*DFThermalDt(r_mjk,t_mjk,beta,a)*DXDisjoint(l,k,pk,mu,beta); //
			else
				res += (2.0*DFThermalDrOnr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta)); //
		}
		if (j!=mk) {
			if (mu==3)
				res += 2.0*DFThermalDt(r_jmk,t_jmk,beta,a)*DXDisjoint(l,j,pj,nu,beta);//
			else
				res += -(2.0*DFThermalDrOnr(r_jmk,t_jmk,beta,a)*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,mk,mu,beta)); //
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
					res += 2.0*DFThermalDt_ij*DXDisjoint(l,i,pi,nu,beta); //
					if (nu==3) 
						res += 2.0*DDFThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDFThermalDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
				}
				else {
					res += 2.0*DFThermalDrOnr_ij*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta); //
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
					res += 2.0*DFThermalDt_ij*DXDisjoint(l,i,pi,mu,beta); //
				}
				else {
					res += - 2.0*DFThermalDrOnr_ij*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,i,j,nu,beta); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==3)
					res += - 2.0*DFThermalDt(r_imj,t_imj,beta,a)*DXDisjoint(l,i,pi,mu,beta);//
				else
					res += 2.0*DFThermalDrOnr(r_imj,t_imj,beta,a)*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,i,k,nu,beta);//
			} 
			if (j==mk && i!=j) {
				if (mu==3)
					res += - 2.0*DFThermalDt(r_ij,t_ij,beta,a)*DXDisjoint(l,i,pi,nu,beta);//
				else
					res += - 2.0*DFThermalDrOnr(r_ij,t_ij,beta,a)*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta);//
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
		res +=  - (4.0*E_jk*DX(l,j,pj,nu)*DX(l,j,k,mu))/pow(a,2) \
				 + (4.0*E_jk*DX(l,j,k,nu)*DX(l,k,pk,mu))/pow(a,2) \
				 - (8.0*E_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(a,4); //
	if (k!=mj)
		res +=  - (4.0*E_mjk*DX(l,mj,k,nu)*DX(l,k,pk,mu))/pow(a,2); //
	if (k!=pj)
		res += + (4.0*E_jmk*DX(l,j,pj,nu)*DX(l,j,mk,mu))/pow(a,2); //
	
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
				res += - (4.0*E_ij*DX(l,i,pi,nu)*DX(l,j,i,mu))/pow(a,2) \
						 - (4.0*E_ij*DX(l,i,pi,mu)*DX(l,j,i,nu))/pow(a,2) \
						 + (8.0*E_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(a,4);
				if (mu==nu)
					res +=  - (4.0*E_ij*T_ij)/pow(a,2); //
			}
			if (k==mj && i!=mj) 
				res +=  (4.0*E_imj*DX(l,i,pi,mu)*DX(l,mj,i,nu))/pow(a,2);
			if (k==pj && i!=j) 
				res += + (4.0*E_ij*DX(l,i,pi,nu)*DX(l,j,i,mu))/pow(a,2);
			
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
		res +=  - (4.0*E_jk*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,k,mu,beta))/pow(a,2) \
				 + (4.0*E_jk*DXDisjoint(l,j,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta))/pow(a,2) \
				 - (8.0*E_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(a,4); //
	if (k!=mj)
		res +=  - (4.0*E_mjk*DXDisjoint(l,mj,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta))/pow(a,2); //
	if (k!=pj)
		res += + (4.0*E_jmk*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,mk,mu,beta))/pow(a,2); //
	
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
				res += - (4.0*E_ij*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta))/pow(a,2) \
						 - (4.0*E_ij*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,j,i,nu,beta))/pow(a,2) \
						 + (8.0*E_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(a,4);
				if (mu==nu)
					res +=  - (4.0*E_ij*T_ij)/pow(a,2); //
			}
			if (k==mj && i!=mj) 
				res +=  (4.0*E_imj*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,mj,i,nu,beta))/pow(a,2);
			if (k==pj && i!=j) 
				res += + (4.0*E_ij*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta))/pow(a,2);
			
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

// mdI_nr<4>
template <> void mdI_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& f, vec& v) {
	if (mu==3) {
		uint nj = negNeigh(j,l.size());
		v[j*4+mu] += -f*DX(l,nj,j,2);
	}
	else if (mu==2) {
		uint pj = posNeigh(j,l.size());
		v[j*4+mu] += -f*DX(l,pj,j,3);
	}
}

// ddI_nr<4>
template <> void ddI_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<4>& l, const number& f, mat& m) {
	if (mu==3 && nu==2) {
		uint nj = negNeigh(j,l.size());
		if (k==nj)
			m(j*4+mu,k*4+nu) += f;
		else if (k==j)
			m(j*4+mu,k*4+nu) -= f;
	}
	else if (mu==2 && nu==3) {
		uint pj = posNeigh(j,l.size());
		if (k==pj)
			m(j*4+mu,k*4+nu) += f;
		else if (k==j)
			m(j*4+mu,k*4+nu) -= f;
	}
}	

// I0<4>
template <> void I0<4> (const uint& j, const Loop<4>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	result += f*(l[j])[2]*DX(l,pj,j,3);
}

// I0Disjoint<4>
template <> void I0Disjoint<4> (const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	result += f*(l[j])[2]*DXDisjoint(l,pj,j,3,beta);
}

// mdIDisjoint_nr
template <> void mdIDisjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& f, vec& v) {
	if (mu==3) {
		uint nj = negNeighDisjoint(j,l.size());
		v[j*4+mu] += -f*DXDisjoint(l,nj,j,2,beta);
	}
	else if (mu==2) {
		uint pj = posNeighDisjoint(j,l.size());
		v[j*4+mu] += -f*DXDisjoint(l,pj,j,3,beta);
	}
}

// ddIDisjoint_nr
template <> void ddIDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& beta, const number& f, mat& m) {
	if (mu==3 && nu==2) {
		uint nj = negNeighDisjoint(j,l.size());
		if (k==nj)
			m(j*4+mu,k*4+nu) += f;
		else if (k==j)
			m(j*4+mu,k*4+nu) -= f;
	}
	else if (mu==2 && nu==3) {
		uint pj = posNeighDisjoint(j,l.size());
		if (k==pj)
			m(j*4+mu,k*4+nu) += f;
		else if (k==j)
			m(j*4+mu,k*4+nu) -= f;
	}
}
		
