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
#include "nrloop2.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		0.0 - static scalar functions
		0.1 - static vector and matrix functions
		1.0 - nr loop scalar functions
		1.1 - nr loop vector and matrix functions
		2 - NEWloopToVector, vectorToLoop
		3 - filename functions
		4 - explicit instatiation, for specific dimensions
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0.0 - static scalar functions
----------------------------------------------------------------------------------------------------------------------------*/

typedef number(*ThermalFunction)(const number&, const number&, const number&, const number&);

typedef number(*VacuumFunction)(const number&, const number&);

// F0
static number F0(const number& x, const number& a) {
	return  -1.0/(4.0*pow(PI,2)*(pow(a,2) + pow(x,2)));
}

// DF0DrOnr
static number DF0DrOnr(const number& x, const number& a) {
	return 1.0/(2.0*pow(PI,2)*pow(pow(a,2) + pow(x,2),2));
}
 
 // DDF0DrDr
static number DDF0DrDr(const number& x, const number& a) {
	return (pow(a,2) - 3.0*pow(x,2))/(2.0*pow(PI,2)*pow(pow(a,2) + pow(x,2),3));
}

// FExp
static number FExp(const number& x, const number& a) {
	return -1.0/(4.0*exp(pow(a,2)/pow(x,2))*pow(PI,2)*pow(x,2));
}

// DFExpDrOnr
static number DFExpDrOnr(const number& x, const number& a) {
	return -(pow(a,2) - pow(x,2))/(2.0*exp(pow(a,2)/pow(x,2))*pow(PI,2)*pow(x,6));
}
 
 // DDFExpDrDr
static number DDFExpDrDr(const number& x, const number& a) {
	return -(pow(a,2) - 3.0*pow(x,2))/(2.0*pow(PI,2)*pow(pow(a,2) + pow(x,2),3));
}

// FDim
static number FDim(const number& x, const number& a) {
	return  - pow(x,-2.0 + a)/(4.0*pow(PI,2));
}

// DFDimDrOnr
static number DFDimDrOnr(const number& x, const number& a) {
	return - ((-2.0 + a)*pow(x,-4.0 + a))/(4.0*pow(PI,2));
}
 
 // DDFDimDrDr
static number DDFDimDrDr(const number& x, const number& a) {
	return  - ((-3.0 + a)*(-2.0 + a)*pow(x,-4.0 + a))/(4.0*pow(PI,2));
}

// G0
static number G0(const number& x, const number& a) {
	return  -1.0/(4.0*exp(pow(x,2)/pow(a,2))*pow(PI,2));
}

// DG0DrOnr
static number DG0DrOnr(const number& x, const number& a) {
	return 1.0/(2.0*pow(a,2)*exp(pow(x,2)/pow(a,2))*pow(PI,2));
}
 
 // DDG0DrDr
static number DDG0DrDr(const number& x, const number& a) {
	return (pow(a,2) - 2.0*pow(x,2))/(2.0*pow(a,4)*exp(pow(x,2)/pow(a,2))*pow(PI,2));
}

// FThermal
static number FThermal(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return  sinh((2.0*PI*ra)/beta)/(4.0*PI*ra*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))); //
}

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
	return (a*sqrt(PI)*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/(exp(pow(r,2)/pow(a,2))*beta)/(-pow(2.0*PI,2));
}

// DGThermalDrOnr
static number DGThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/(a*exp(pow(r,2)/pow(a,2))*beta)/(-pow(2.0*PI,2));
}

// DGThermalDt
static number DGThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	return (a*sqrt(PI)*jbderivutheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (exp(pow(r,2)/pow(a,2))*pow(beta,2))/(-pow(2.0*PI,2));
}

// DDGThermalDrDr
static number DDGThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*(pow(a,2) - 2.0*pow(r,2))*jbdtheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (pow(a,3)*exp(pow(r,2)/pow(a,2))*beta)/(-pow(2.0*PI,2));
}

// DDGThermalDtDrOnr
static number DDGThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (-2.0*sqrt(PI)*jbderivutheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (a*exp(pow(r,2)/pow(a,2))*pow(beta,2))/(-pow(2.0*PI,2));
}

// DDGThermalDtDt
static number DDGThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	return (a*sqrt(PI)*jbderiv2utheta3(t/beta,exp(-(pow(a,2)*pow(PI,2))/pow(beta,2))))/
   (exp(pow(r,2)/pow(a,2))*pow(beta,3))/(-pow(2.0*PI,2));
}

// G2Thermal
static number G2Thermal(const number& r, const number& t, const number& beta, const number& a) {
	return  exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))/(-pow(2.0*PI,2));
}

// DG2ThermalDrOnr
static number DG2ThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return ((-2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2))))/pow(a,2))/(-pow(2.0*PI,2));
}

// DG2ThermalDt
static number DG2ThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	return - ((exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*beta*sin((2.0*PI*t)/beta))/\
 (pow(a,2)*PI))/(-pow(2.0*PI,2));
}

// DDG2ThermalDrDr
static number DDG2ThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	return ((-2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2))))/pow(a,2) \
 + (4.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*pow(r,2))/pow(a,4))/(-pow(2.0*PI,2));
}

// DDG2ThermalDtDrOnr
static number DDG2ThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return  ((2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*beta*sin((2.0*PI*t)/beta))/\
 (pow(a,4)*PI))/(-pow(2.0*PI,2));
}

// DDG2ThermalDtDt
static number DDG2ThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	return  ((-2.0*exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*cos((2.0*PI*t)/beta))/\
 pow(a,2) + (exp(-(pow(r,2)/pow(a,2)) + (pow(beta,2)*(-1 + cos((2.0*PI*t)/beta)))/(2.0*pow(a,2)*pow(PI,2)))*pow(beta,2)*\
 pow(sin((2.0*PI*t)/beta),2))/(pow(a,4)*pow(PI,2)))/(-pow(2.0*PI,2));
}

// V0Generic
template <uint Dim>
static void V0Generic (const uint& j, const uint& k, const Loop<Dim>& l, VacuumFunction FVacuum, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = posNeigh(j,l.size());
		uint pk = posNeigh(k,l.size());
		number x_jk = Distance(l[j],l[k]);
		
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*FVacuum(x_jk,a);	
	}
	
}

// V0DisjointGeneric
template <uint Dim>
static void V0DisjointGeneric (const uint& j, const uint& k, const Loop<Dim>& l, VacuumFunction FVacuum, const number& beta, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = posNeighDisjoint(j,l.size());
		uint pk = posNeighDisjoint(k,l.size());
		number x_jk = DistanceDisjoint(l[j],l[k],beta);
	
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*DotDisjoint(l[pj],l[j],l[pk],l[k],beta)*FVacuum(x_jk,a);	
	}
	
}

// V0DisjointLRGeneric
template <uint Dim>
static void V0DisjointLRGeneric (const uint& j, const uint& k, const Loop<Dim>& l, VacuumFunction FVacuum, const number& beta, const number& a, const number& f, number& result) {

	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		V0DisjointGeneric(j,k,l,FVacuum,beta,a,f,result);
	
}

// VthGeneric
template <uint Dim>
static void VthGeneric (const uint& j, const uint& k, const Loop<Dim>& l, ThermalFunction Fth, const number& beta, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = posNeigh(j,l.size());
		uint pk = posNeigh(k,l.size());
		number r = SpatialDistance(l[j],l[k]);
		number t = DX(l,k,j,Dim-1);
		
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*Fth(r,t,beta,a);	
	}
	
}

// VthDisjointGeneric
template <uint Dim>
static void VthDisjointGeneric (const uint& j, const uint& k, const Loop<Dim>& l, ThermalFunction Fth, const number& beta, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = posNeighDisjoint(j,l.size());
		uint pk = posNeighDisjoint(k,l.size());
		
		number r = SpatialDistance(l[j],l[k]);
		number t = DXDisjoint(l,k,j,Dim-1,beta);
	
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*DotDisjoint(l[pj],l[j],l[pk],l[k],beta)*Fth(r,t,beta,a);	
	}
	
}

// VthDisjointLRGeneric
template <uint Dim>
static void VthDisjointLRGeneric (const uint& j, const uint& k, const Loop<Dim>& l, ThermalFunction Fth, const number& beta, const number& a, const number& f, number& result) {

	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		VthDisjointGeneric(j,k,l,Fth,beta,a,f,result);
	
}

// VthDisjointLR2Generic
template <uint Dim>
static void VthDisjointLR2Generic (const uint& j, const uint& k, const Loop<Dim>& l, ThermalFunction Fth, const number& beta, const number& a, const number& f, number& result) {

	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		VthDisjointGeneric(j,k,l,Fth,beta,a,f,result);
	else
		VthDisjointGeneric(j,k,l,Fth,beta,0.0,f,result);
	
}

/*----------------------------------------------------------------------------------------------------------------------------
	0.1 - static vector and matrix functions
----------------------------------------------------------------------------------------------------------------------------*/
 
// mdV0Generic
template<uint Dim>
static void mdV0Generic(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					VacuumFunction FVacuum, VacuumFunction dFVacuumdronr, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number x_ij = Distance(l[i],l[j]);
		number F0_ij = FVacuum(x_ij,a);
		number DF0DrOnr_ij = dFVacuumdronr(x_ij,a);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		
		res += 2.0*F0_ij*DX(l,i,pi,mu) \
 				+ 2.0*DF0DrOnr_ij*DX(l,j,i,mu)*T_ij;
	}
	
	if (i!=mj) {
		number x_imj = Distance(l[i],l[mj]);
		number F0_imj = FVacuum(x_imj,a);
		res +=  - 2.0*F0_imj*DX(l,i,pi,mu); //
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (i==j)
		res += 2.0*(l[j])[mu]*FVacuum(0.0,a); 
	if (i==mj)
		res += (-(l[mj])[mu]*FVacuum(0.0,a));
	if (i==pj)
		res += (-(l[pj])[mu]*FVacuum(0.0,a));
		
	v[j*Dim+mu] += -f*(-pow(2.0*PI,2))*res;
}

// ddV0Generic
template<uint Dim>
static void ddV0Generic(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						VacuumFunction FVacuum, VacuumFunction dFVacuumdronr,VacuumFunction ddFVacuumdrdr,\
						 const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeigh(j,l.size());
	uint pj = posNeigh(j,l.size());	
	uint mk = negNeigh(k,l.size());
	uint pk = posNeigh(k,l.size());
	
	number x_jk = Distance(l[j],l[k]);
	number x_jmk = Distance(l[j],l[mk]);
	number x_mjk = Distance(l[mj],l[k]);
	number x_mjmk = Distance(l[mj],l[mk]);
	
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*FVacuum(x_mjmk,a) \
			 + 2.0*FVacuum(x_jk,a) \
			 - 2.0*dFVacuumdronr(x_jk,a)*T_jk;
		if (k!=mj)
			res += - 2.0*FVacuum(x_mjk,a);
		if (k!=pj)
			res += - 2.0*FVacuum(x_jmk,a);
	}

	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res +=  + 2.0*dFVacuumdronr(x_jk,a)*DX(l,j,pj,nu)*DX(l,j,k,mu) \
				 - 2.0*dFVacuumdronr(x_jk,a)*DX(l,j,k,nu)*DX(l,k,pk,mu) \
				 - (2.0*ddFVacuumdrdr(x_jk,a)*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(x_jk,2) \
				 + (2.0*dFVacuumdronr(x_jk,a)*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(x_jk,2);
	if (k!=mj)
		res += + 2.0*dFVacuumdronr(x_mjk,a)*DX(l,mj,k,nu)*DX(l,k,pk,mu);
	if (k!=pj)
		res += - 2.0*dFVacuumdronr(x_jmk,a)*DX(l,j,pj,nu)*DX(l,j,mk,mu);
		
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number x_ij, x_imj, T_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = (i==(l.size()-1)? 0: i+1);
			x_ij = Distance(l[i],l[j]);
			x_imj = Distance(l[i],l[mj]);
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);
			
			if (k==j && i!=j) {
				res +=  2.0*dFVacuumdronr(x_ij,a)*DX(l,i,pi,mu)*DX(l,j,i,nu) \
						 + 2.0*dFVacuumdronr(x_ij,a)*DX(l,i,pi,nu)*DX(l,j,i,mu) \
						 + (2.0*ddFVacuumdrdr(x_ij,a)*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(x_ij,2) \
						 - (2.0*dFVacuumdronr(x_ij,a)*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(x_ij,2);
				if (mu==nu)
					res += + 2.0*dFVacuumdronr(x_ij,a)*T_ij; //
			}
			if (k==mj && i!=mj) 
				res += 2.0*dFVacuumdronr(x_imj,a)*DX(l,i,pi,mu)*DX(l,i,k,nu);
			if (k==pj && i!=j) 
				res += - 2.0*dFVacuumdronr(x_ij,a)*DX(l,i,pi,nu)*DX(l,j,i,mu);
			
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += (2.0*FVacuum(0.0,a));
	if (k==mj && mu==nu)
		res += (-1.0*FVacuum(0.0,a));
	if (k==pj && mu==nu)
		res += (-1.0*FVacuum(0.0,a));
	
	m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
	
}

// mdVthGeneric
template<uint Dim>
static void mdVthGeneric(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, \
					 	const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DX(l,j,i,Dim-1); // checked order
		number FTH_ij = Fth(r_ij,t_ij,beta,a);
		number DFTHDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
		number DFTHDt_ij = dFthdt(r_ij,t_ij,beta,a);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		
		res += 2.0*FTH_ij*(-DX(l,i,mu));
		if (mu<(Dim-1))
			res += 2.0*DFTHDrOnr_ij*DX(l,j,i,mu)*T_ij;
		else
			res += 2.0*DFTHDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DX(l,mj,i,Dim-1); // checked order
		number FTH_imj = Fth(r_imj,t_imj,beta,a);
		res +=  -2.0*FTH_imj*(-DX(l,i,mu)); //
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (i==j)
		res += 2.0*(l[j])[mu]*Fth(0.0,0.0,beta,a); 
	if (i==mj)
		res += (-(l[mj])[mu]*Fth(0.0,0.0,beta,a));
	if (i==pj)
		res += (-(l[pj])[mu]*Fth(0.0,0.0,beta,a));
		
	v[j*Dim+mu] += -f*(-pow(2.0*PI,2))*res;
}

// ddVthGeneric
template <uint Dim>
static void ddVthGeneric(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
	ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, ThermalFunction ddFthdrdr, ThermalFunction ddFthdrdtonr, ThermalFunction ddFthdtdt,\
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
	
	number t_jk = DX(l,k,j,Dim-1);
	number t_mjk = DX(l,k,mj,Dim-1);
	number t_jmk = DX(l,mk,j,Dim-1);
	number t_mjmk = DX(l,mk,mj,Dim-1);
		
	number T_jk = Dot(l[pj],l[j],l[pk],l[k]);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*Fth(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*Fth(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*Fth(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*Fth(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DFTHDrOnr_jk = dFthdronr(r_jk,t_jk,beta,a);
			number DFTHDt_jk = dFthdt(r_jk,t_jk,beta,a);
			number DDFTHDrDr_jk = ddFthdrdr(r_jk,t_jk,beta,a);
			number DDFTHDrDtOnr_jk = ddFthdrdtonr(r_jk,t_jk,beta,a);
			number DDFTHDtDt_jk = ddFthdtdt(r_jk,t_jk,beta,a);
			if (mu==(Dim-1)) {
				res += - 2.0*DFTHDt_jk*(-DX(l,j,nu)); //
				if (nu<(Dim-1))
					res += (2.0*DDFTHDrDtOnr_jk*DX(l,j,k,nu)*T_jk); //
				else
					res += -2.0*DDFTHDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFTHDrOnr_jk*(-DX(l,j,nu))*DX(l,j,k,mu); //
				if (nu==(Dim-1)) 
					res += + 2.0*DDFTHDrDtOnr_jk*DX(l,j,k,mu)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFTHDrOnr+DDFTHDrDr->0 as r->0
					res += (2.0*DFTHDrOnr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFTHDrDr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFTHDrOnr_jk*T_jk);//
			}
			
			if (nu==(Dim-1)) {
				res += 2.0*DFTHDt_jk*(-DX(l,k,mu)); //
			}
			else {
				res += - 2.0*DFTHDrOnr_jk*DX(l,j,k,nu)*(-DX(l,k,mu)); //
			}
			
		}
		if (k!=mj) {
			if (nu==(Dim-1))
				res += - 2.0*dFthdt(r_mjk,t_mjk,beta,a)*(-DX(l,k,mu)); //
			else
				res += 2.0*dFthdronr(r_mjk,t_mjk,beta,a)*DX(l,mj,k,nu)*(-DX(l,k,mu)); //
		}
		if (j!=mk) {
			if (mu==(Dim-1))
				res += 2.0*dFthdt(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu));//
			else
				res += -(2.0*dFthdronr(r_jmk,t_jmk,beta,a)*(-DX(l,j,nu))*DX(l,j,mk,mu)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DFTHDrOnr_ij;	
		number DFTHDt_ij;	
		number DDFTHDrDr_ij;	
		number DDFTHDrDtOnr_ij;
		number DDFTHDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeigh(i,l.size());
			T_ij = Dot(l[pi],l[i],l[pj],l[j]);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DX(l,j,i,Dim-1);
			t_imj = DX(l,mj,i,Dim-1);

			if (k==j && i!=j) {
				DFTHDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
				DFTHDt_ij = dFthdt(r_ij,t_ij,beta,a);
				DDFTHDrDr_ij = ddFthdrdr(r_ij,t_ij,beta,a);
				DDFTHDrDtOnr_ij = ddFthdrdtonr(r_ij,t_ij,beta,a);
				DDFTHDtDt_ij = ddFthdtdt(r_ij,t_ij,beta,a);
				if (mu==(Dim-1)) {
					res += 2.0*DFTHDt_ij*(-DX(l,i,nu)); //
					if (nu==(Dim-1)) 
						res += 2.0*DDFTHDtDt_ij*T_ij; //
					else
						res += 2.0*DDFTHDrDtOnr_ij*DX(l,j,i,nu)*T_ij; //
				}
				else {
					res += 2.0*DFTHDrOnr_ij*(-DX(l,i,nu))*DX(l,j,i,mu); //
					if (nu==(Dim-1))
						res += 2.0*DDFTHDrDtOnr_ij*DX(l,j,i,mu)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DFTHDrOnr+DDFTHDrDr->0 as r->0
						res += - (2.0*DFTHDrOnr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDFTHDrDr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DFTHDrOnr_ij*T_ij; //
				}
				if (nu==(Dim-1)) {
					res += 2.0*DFTHDt_ij*(-DX(l,i,mu)); //
				}
				else {
					res += - 2.0*DFTHDrOnr_ij*(-DX(l,i,mu))*DX(l,i,j,nu); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==(Dim-1))
					res += - 2.0*dFthdt(r_imj,t_imj,beta,a)*(-DX(l,i,mu));//
				else
					res += 2.0*dFthdronr(r_imj,t_imj,beta,a)*(-DX(l,i,mu))*DX(l,i,k,nu);//
			} 
			if (j==mk && i!=j) {
				if (mu==(Dim-1))
					res += - 2.0*dFthdt(r_ij,t_ij,beta,a)*(-DX(l,i,nu));//
				else
					res += - 2.0*dFthdronr(r_ij,t_ij,beta,a)*(-DX(l,i,nu))*DX(l,j,i,mu);//
			}	
		}		
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (k==j && mu==nu)
		res += 2.0*Fth(0.0,0.0,beta,a);
	if (k==mj && mu==nu)
		res += (-1.0*Fth(0.0,0.0,beta,a));
	if (k==pj && mu==nu)
		res += (-1.0*Fth(0.0,0.0,beta,a));
	
	m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
}

// mdVthDisjointGeneric
template<uint Dim>
static void mdVthDisjointGeneric(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, \
					 	const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeighDisjoint(j,l.size());
	uint mj = negNeighDisjoint(j,l.size());
	uint pi = posNeighDisjoint(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DXDisjoint(l,j,i,(Dim-1),beta); // checked order
		number FTH_ij = Fth(r_ij,t_ij,beta,a);
		number DFTHDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
		number DFTHDt_ij = dFthdt(r_ij,t_ij,beta,a);
		number T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
		
		res += 2.0*FTH_ij*(-DXDisjoint(l,i,mu,beta));
		if (mu<(Dim-1))
			res += 2.0*DFTHDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij;
		else
			res += 2.0*DFTHDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DXDisjoint(l,mj,i,(Dim-1),beta); // checked order
		number FTH_imj = Fth(r_imj,t_imj,beta,a);
		res +=  -2.0*FTH_imj*(-DXDisjoint(l,i,mu,beta)); //
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (mu==(Dim-1)) {
		if (i==j) { // bit of a fudge in terms of sum over i
			res += (DXDisjoint(l,j,mj,mu,beta) + DXDisjoint(l,j,pj,mu,beta))*Fth(0.0,0.0,beta,a);
		}
	}
	else {
		if (i==j)
			res += 2.0*(l[j])[mu]*Fth(0.0,0.0,beta,a);
		if (i==mj)
			res += (-(l[mj])[mu]*Fth(0.0,0.0,beta,a));
		if (i==pj)
			res += (-(l[pj])[mu]*Fth(0.0,0.0,beta,a));
	}
		
	v[j*Dim+mu] += -f*(-pow(2.0*PI,2))*res;
}

// ddVthDisjointGeneric
template <uint Dim>
static void ddVthDisjointGeneric(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
	ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, ThermalFunction ddFthdrdr, ThermalFunction ddFthdrdtonr, ThermalFunction ddFthdtdt,\
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
	
	number t_jk = DXDisjoint(l,k,j,(Dim-1),beta);
	number t_mjk = DXDisjoint(l,k,mj,(Dim-1),beta);
	number t_jmk = DXDisjoint(l,mk,j,(Dim-1),beta);
	number t_mjmk = DXDisjoint(l,mk,mj,(Dim-1),beta);
		
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*Fth(r_mjmk,t_mjmk,beta,a)\
					+ 2.0*Fth(r_jk,t_jk,beta,a);//
		if (k!=mj)
			res +=  - 2.0*Fth(r_mjk,t_mjk,beta,a);//
		if (k!=pj)
			res += - 2.0*Fth(r_jmk,t_jmk,beta,a);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DFTHDrOnr_jk = dFthdronr(r_jk,t_jk,beta,a);
			number DFTHDt_jk = dFthdt(r_jk,t_jk,beta,a);
			number DDFTHDrDr_jk = ddFthdrdr(r_jk,t_jk,beta,a);
			number DDFTHDrDtOnr_jk = ddFthdrdtonr(r_jk,t_jk,beta,a);
			number DDFTHDtDt_jk = ddFthdtdt(r_jk,t_jk,beta,a);
			if (mu==(Dim-1)) {
				res += - 2.0*DFTHDt_jk*(-DXDisjoint(l,j,nu,beta)); //
				if (nu<(Dim-1))
					res += (2.0*DDFTHDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
				else
					res += -2.0*DDFTHDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFTHDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
				if (nu==(Dim-1)) 
					res += + 2.0*DDFTHDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFTHDrOnr+DDFTHDrDr->0 as r->0
					res += (2.0*DFTHDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFTHDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFTHDrOnr_jk*T_jk);//
			}
			
			if (nu==(Dim-1)) {
				res += 2.0*DFTHDt_jk*(-DXDisjoint(l,k,mu,beta)); //
			}
			else {
				res += - 2.0*DFTHDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			
		}
		if (k!=mj) {
			if (nu==(Dim-1))
				res += - 2.0*dFthdt(r_mjk,t_mjk,beta,a)*(-DXDisjoint(l,k,mu,beta)); //
			else
				res += 2.0*dFthdronr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
		}
		if (j!=mk) {
			if (mu==(Dim-1))
				res += 2.0*dFthdt(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta));//
			else
				res += -(2.0*dFthdronr(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DFTHDrOnr_ij;	
		number DFTHDt_ij;	
		number DDFTHDrDr_ij;	
		number DDFTHDrDtOnr_ij;
		number DDFTHDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeighDisjoint(i,l.size());
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DXDisjoint(l,j,i,(Dim-1),beta);
			t_imj = DXDisjoint(l,mj,i,(Dim-1),beta);

			if (k==j && i!=j) {
				DFTHDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
				DFTHDt_ij = dFthdt(r_ij,t_ij,beta,a);
				DDFTHDrDr_ij = ddFthdrdr(r_ij,t_ij,beta,a);
				DDFTHDrDtOnr_ij = ddFthdrdtonr(r_ij,t_ij,beta,a);
				DDFTHDtDt_ij = ddFthdtdt(r_ij,t_ij,beta,a);
				if (mu==(Dim-1)) {
					res += 2.0*DFTHDt_ij*(-DXDisjoint(l,i,nu,beta)); //
					if (nu==(Dim-1)) 
						res += 2.0*DDFTHDtDt_ij*T_ij; //
					else
						res += 2.0*DDFTHDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
				}
				else {
					res += 2.0*DFTHDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
					if (nu==(Dim-1))
						res += 2.0*DDFTHDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DFTHDrOnr+DDFTHDrDr->0 as r->0
						res += - (2.0*DFTHDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDFTHDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DFTHDrOnr_ij*T_ij; //
				}
				if (nu==(Dim-1)) {
					res += 2.0*DFTHDt_ij*(-DXDisjoint(l,i,mu,beta)); //
				}
				else {
					res += - 2.0*DFTHDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==(Dim-1))
					res += - 2.0*dFthdt(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta));//
				else
					res += 2.0*dFthdronr(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
			} 
			if (j==mk && i!=j) {
				if (mu==(Dim-1))
					res += - 2.0*dFthdt(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta));//
				else
					res += - 2.0*dFthdronr(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
			}	
		}		
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (k==j && mu==nu)
		res += 2.0*Fth(0.0,0.0,beta,a);
	if (k==mj && mu==nu)
		res += (-1.0*Fth(0.0,0.0,beta,a));
	if (k==pj && mu==nu)
		res += (-1.0*Fth(0.0,0.0,beta,a));
	
	m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
}

// mdVthDisjointLRGeneric
template<uint Dim>
static void mdVthDisjointLRGeneric(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && i<N/2) || (j>=N/2 && i>=N/2))
		mdVthDisjointGeneric(j,mu,i,l,Fth,dFthdronr,dFthdt,beta,a,f,v);
}

// ddVthDisjointLRGeneric
template <uint Dim>
static void ddVthDisjointLRGeneric(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
	ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, ThermalFunction ddFthdrdr, ThermalFunction ddFthdrdtonr, ThermalFunction ddFthdtdt,\
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
	
		number t_jk = DXDisjoint(l,k,j,(Dim-1),beta);
		number t_mjk = DXDisjoint(l,k,mj,(Dim-1),beta);
		number t_jmk = DXDisjoint(l,mk,j,(Dim-1),beta);
		number t_mjmk = DXDisjoint(l,mk,mj,(Dim-1),beta);
		
		number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
		// terms where mu==nu, without sums
		if (mu==nu) {
			if (k!=j)
				res += 2.0*Fth(r_mjmk,t_mjmk,beta,a)\
						+ 2.0*Fth(r_jk,t_jk,beta,a);//
			if (k!=mj)
				res +=  - 2.0*Fth(r_mjk,t_mjk,beta,a);//
			if (k!=pj)
				res += - 2.0*Fth(r_jmk,t_jmk,beta,a);//
		}

		// terms where mu not nexcessarily equal to nu, without sums
			if (k!=j) {
				number DFTHDrOnr_jk = dFthdronr(r_jk,t_jk,beta,a);
				number DFTHDt_jk = dFthdt(r_jk,t_jk,beta,a);
				number DDFTHDrDr_jk = ddFthdrdr(r_jk,t_jk,beta,a);
				number DDFTHDrDtOnr_jk = ddFthdrdtonr(r_jk,t_jk,beta,a);
				number DDFTHDtDt_jk = ddFthdtdt(r_jk,t_jk,beta,a);
				if (mu==(Dim-1)) {
					res += - 2.0*DFTHDt_jk*(-DXDisjoint(l,j,nu,beta)); //
					if (nu<(Dim-1))
						res += (2.0*DDFTHDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
					else
						res += -2.0*DDFTHDtDt_jk*T_jk;//
				}
				else {
					res += + 2.0*DFTHDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
					if (nu==(Dim-1)) 
						res += + 2.0*DDFTHDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
					else if (r_jk>MIN_NUMBER) //as -DFTHDrOnr+DDFTHDrDr->0 as r->0
						res += (2.0*DFTHDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
								- (2.0*DDFTHDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
					if (nu==mu)
						res += - (2.0*DFTHDrOnr_jk*T_jk);//
				}
			
				if (nu==(Dim-1)) {
					res += 2.0*DFTHDt_jk*(-DXDisjoint(l,k,mu,beta)); //
				}
				else {
					res += - 2.0*DFTHDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
				}
			
			}
			if (k!=mj) {
				if (nu==(Dim-1))
					res += - 2.0*dFthdt(r_mjk,t_mjk,beta,a)*(-DXDisjoint(l,k,mu,beta)); //
				else
					res += 2.0*dFthdronr(r_mjk,t_mjk,beta,a)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			if (j!=mk) {
				if (mu==(Dim-1))
					res += 2.0*dFthdt(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta));//
				else
					res += -(2.0*dFthdronr(r_jmk,t_jmk,beta,a)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
			}
	
		// terms with sums
		if (k==j || k==mj || k==pj) {
	
			uint pi;
			number T_ij;
			number r_ij, r_imj;
			number t_ij, t_imj;
			number DFTHDrOnr_ij;	
			number DFTHDt_ij;	
			number DDFTHDrDr_ij;	
			number DDFTHDrDtOnr_ij;
			number DDFTHDtDt_ij;
			
			// for LR
			uint iMin = (j<N/2? 0: N/2);
			uint iMax = (j<N/2? N/2: N);
		
			for (uint i=iMin; i<iMax; i++) {
		
				pi = posNeighDisjoint(i,l.size());
				T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
				r_ij = SpatialDistance(l[j],l[i]);
				r_imj = SpatialDistance(l[i],l[mj]);	
				t_ij = DXDisjoint(l,j,i,(Dim-1),beta);
				t_imj = DXDisjoint(l,mj,i,(Dim-1),beta);

				if (k==j && i!=j) {
					DFTHDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
					DFTHDt_ij = dFthdt(r_ij,t_ij,beta,a);
					DDFTHDrDr_ij = ddFthdrdr(r_ij,t_ij,beta,a);
					DDFTHDrDtOnr_ij = ddFthdrdtonr(r_ij,t_ij,beta,a);
					DDFTHDtDt_ij = ddFthdtdt(r_ij,t_ij,beta,a);
					if (mu==(Dim-1)) {
						res += 2.0*DFTHDt_ij*(-DXDisjoint(l,i,nu,beta)); //
						if (nu==(Dim-1)) 
							res += 2.0*DDFTHDtDt_ij*T_ij; //
						else
							res += 2.0*DDFTHDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
					}
					else {
						res += 2.0*DFTHDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
						if (nu==(Dim-1))
							res += 2.0*DDFTHDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
						else if (r_ij>MIN_NUMBER) { // as -DFTHDrOnr+DDFTHDrDr->0 as r->0
							res += - (2.0*DFTHDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
								+ (2.0*DDFTHDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
						}
						if (nu==mu)
							res += 2.0*DFTHDrOnr_ij*T_ij; //
					}
					if (nu==(Dim-1)) {
						res += 2.0*DFTHDt_ij*(-DXDisjoint(l,i,mu,beta)); //
					}
					else {
						res += - 2.0*DFTHDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
					}
				} 
				if (j==pk && i!=mj) {
					if (nu==(Dim-1))
						res += - 2.0*dFthdt(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta));//
					else
						res += 2.0*dFthdronr(r_imj,t_imj,beta,a)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
				} 
				if (j==mk && i!=j) {
					if (mu==(Dim-1))
						res += - 2.0*dFthdt(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta));//
					else
						res += - 2.0*dFthdronr(r_ij,t_ij,beta,a)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
				}	
			}		
		}
	
		//coincident terms
		// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
		if (k==j && mu==nu)
			res += 2.0*Fth(0.0,0.0,beta,a);
		if (k==mj && mu==nu)
			res += (-1.0*Fth(0.0,0.0,beta,a));
		if (k==pj && mu==nu)
			res += (-1.0*Fth(0.0,0.0,beta,a));
	
		m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
	}
}

// mdVthDisjointLR2Generic
template<uint Dim>
static void mdVthDisjointLR2Generic(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, \
					 	const number& beta, const number& a, const number& f, vec& v)  {
	uint N = l.size();
	if ((j<N/2 && i<N/2) || (j>=N/2 && i>=N/2))
		mdVthDisjointGeneric(j,mu,i,l,Fth,dFthdronr,dFthdt,beta,a,f,v);
	else
		mdVthDisjointGeneric(j,mu,i,l,Fth,dFthdronr,dFthdt,beta,0.0,f,v);
}

// ddVthDisjointLR2Generic
template<uint Dim>
static void ddVthDisjointLR2Generic(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
	ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, ThermalFunction ddFthdrdr, ThermalFunction ddFthdrdtonr, ThermalFunction ddFthdtdt,\
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
	
	number t_jk = DXDisjoint(l,k,j,(Dim-1),beta);
	number t_mjk = DXDisjoint(l,k,mj,(Dim-1),beta);
	number t_jmk = DXDisjoint(l,mk,j,(Dim-1),beta);
	number t_mjmk = DXDisjoint(l,mk,mj,(Dim-1),beta);
		
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*Fth(r_mjmk,t_mjmk,beta,a_jk)\
					+ 2.0*Fth(r_jk,t_jk,beta,a_jk);//
		if (k!=mj)
			res +=  - 2.0*Fth(r_mjk,t_mjk,beta,a_jk);//
		if (k!=pj)
			res += - 2.0*Fth(r_jmk,t_jmk,beta,a_jk);//
	}

	// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j) {
			number DFTHDrOnr_jk = dFthdronr(r_jk,t_jk,beta,a_jk);
			number DFTHDt_jk = dFthdt(r_jk,t_jk,beta,a_jk);
			number DDFTHDrDr_jk = ddFthdrdr(r_jk,t_jk,beta,a_jk);
			number DDFTHDrDtOnr_jk = ddFthdrdtonr(r_jk,t_jk,beta,a_jk);
			number DDFTHDtDt_jk = ddFthdtdt(r_jk,t_jk,beta,a_jk);
			if (mu==(Dim-1)) {
				res += - 2.0*DFTHDt_jk*(-DXDisjoint(l,j,nu,beta)); //
				if (nu<(Dim-1))
					res += (2.0*DDFTHDrDtOnr_jk*DXDisjoint(l,j,k,nu,beta)*T_jk); //
				else
					res += -2.0*DDFTHDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFTHDrOnr_jk*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,k,mu,beta); //
				if (nu==(Dim-1)) 
					res += + 2.0*DDFTHDrDtOnr_jk*DXDisjoint(l,j,k,mu,beta)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFTHDrOnr+DDFTHDrDr->0 as r->0
					res += (2.0*DFTHDrOnr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFTHDrDr_jk*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFTHDrOnr_jk*T_jk);//
			}
			
			if (nu==(Dim-1)) {
				res += 2.0*DFTHDt_jk*(-DXDisjoint(l,k,mu,beta)); //
			}
			else {
				res += - 2.0*DFTHDrOnr_jk*DXDisjoint(l,j,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
			}
			
		}
		if (k!=mj) {
			if (nu==(Dim-1))
				res += - 2.0*dFthdt(r_mjk,t_mjk,beta,a_jk)*(-DXDisjoint(l,k,mu,beta)); //
			else
				res += 2.0*dFthdronr(r_mjk,t_mjk,beta,a_jk)*DXDisjoint(l,mj,k,nu,beta)*(-DXDisjoint(l,k,mu,beta)); //
		}
		if (j!=mk) {
			if (mu==(Dim-1))
				res += 2.0*dFthdt(r_jmk,t_jmk,beta,a_jk)*(-DXDisjoint(l,j,nu,beta));//
			else
				res += -(2.0*dFthdronr(r_jmk,t_jmk,beta,a_jk)*(-DXDisjoint(l,j,nu,beta))*DXDisjoint(l,j,mk,mu,beta)); //
		}
	
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
			
		number a_ij;
		uint pi;
		number T_ij;
		number r_ij, r_imj;
		number t_ij, t_imj;
		number DFTHDrOnr_ij;	
		number DFTHDt_ij;	
		number DDFTHDrDr_ij;	
		number DDFTHDrDtOnr_ij;
		number DDFTHDtDt_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = posNeighDisjoint(i,l.size());
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);			
			r_ij = SpatialDistance(l[j],l[i]);
			r_imj = SpatialDistance(l[i],l[mj]);	
			t_ij = DXDisjoint(l,j,i,(Dim-1),beta);
			t_imj = DXDisjoint(l,mj,i,(Dim-1),beta);
			a_ij = (((j<N/2 && i<N/2) || (j>=N/2 && i>=N/2))? a:0.0);

			if (k==j && i!=j) {
				DFTHDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a_ij);
				DFTHDt_ij = dFthdt(r_ij,t_ij,beta,a_ij);
				DDFTHDrDr_ij = ddFthdrdr(r_ij,t_ij,beta,a_ij);
				DDFTHDrDtOnr_ij = ddFthdrdtonr(r_ij,t_ij,beta,a_ij);
				DDFTHDtDt_ij = ddFthdtdt(r_ij,t_ij,beta,a_ij);
				if (mu==(Dim-1)) {
					res += 2.0*DFTHDt_ij*(-DXDisjoint(l,i,nu,beta)); //
					if (nu==(Dim-1)) 
						res += 2.0*DDFTHDtDt_ij*T_ij; //
					else
						res += 2.0*DDFTHDrDtOnr_ij*DXDisjoint(l,j,i,nu,beta)*T_ij; //
				}
				else {
					res += 2.0*DFTHDrOnr_ij*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta); //
					if (nu==(Dim-1))
						res += 2.0*DDFTHDrDtOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DFTHDrOnr+DDFTHDrDr->0 as r->0
						res += - (2.0*DFTHDrOnr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDFTHDrDr_ij*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DFTHDrOnr_ij*T_ij; //
				}
				if (nu==(Dim-1)) {
					res += 2.0*DFTHDt_ij*(-DXDisjoint(l,i,mu,beta)); //
				}
				else {
					res += - 2.0*DFTHDrOnr_ij*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,j,nu,beta); //
				}
			} 
			if (j==pk && i!=mj) {
				if (nu==(Dim-1))
					res += - 2.0*dFthdt(r_imj,t_imj,beta,a_ij)*(-DXDisjoint(l,i,mu,beta));//
				else
					res += 2.0*dFthdronr(r_imj,t_imj,beta,a_ij)*(-DXDisjoint(l,i,mu,beta))*DXDisjoint(l,i,k,nu,beta);//
			} 
			if (j==mk && i!=j) {
				if (mu==(Dim-1))
					res += - 2.0*dFthdt(r_ij,t_ij,beta,a_ij)*(-DXDisjoint(l,i,nu,beta));//
				else
					res += - 2.0*dFthdronr(r_ij,t_ij,beta,a_ij)*(-DXDisjoint(l,i,nu,beta))*DXDisjoint(l,j,i,mu,beta);//
			}	
		}		
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (k==j && mu==nu)
		res += 2.0*Fth(0.0,0.0,beta,a);
	if (k==mj && mu==nu)
		res += (-1.0*Fth(0.0,0.0,beta,a));
	if (k==pj && mu==nu)
		res += (-1.0*Fth(0.0,0.0,beta,a));
	
	m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
}

// mdV0DisjointGeneric
template<uint Dim>
static void mdV0DisjointGeneric(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					VacuumFunction FVacuum, VacuumFunction dFVacuumdronr, const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeighDisjoint(j,l.size());
	uint mj = negNeighDisjoint(j,l.size());
	uint pi = posNeighDisjoint(i,l.size());
		
	if (i!=j) {
		number x_ij = DistanceDisjoint(l[i],l[j],beta);
		number F0_ij = FVacuum(x_ij,a);
		number DF0DrOnr_ij = dFVacuumdronr(x_ij,a);
		number T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
		
		res += 2.0*F0_ij*DXDisjoint(l,i,pi,mu,beta) \
 				+ 2.0*DF0DrOnr_ij*DXDisjoint(l,j,i,mu,beta)*T_ij;
	}
	
	if (i!=mj) {
		number x_imj = DistanceDisjoint(l[i],l[mj],beta);
		number F0_imj = FVacuum(x_imj,a);
		res +=  - 2.0*F0_imj*DXDisjoint(l,i,pi,mu,beta); //
	}
	
	//coincident terms
	// extra factor of (-1.0/pow(2.0*PI,2)) due to the fact that we are treating the green's function here
	if (i==j)
		res += 2.0*(l[j])[mu]*FVacuum(0.0,a); 
	if (i==mj)
		res += (-(l[mj])[mu]*FVacuum(0.0,a));
	if (i==pj)
		res += (-(l[pj])[mu]*FVacuum(0.0,a));
		
	v[j*Dim+mu] += -f*(-pow(2.0*PI,2))*res;

}

// ddV0DisjointGeneric
template<uint Dim>
static void ddV0DisjointGeneric(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						VacuumFunction FVacuum, VacuumFunction dFVacuumdronr,VacuumFunction ddFVacuumdrdr,\
						 const number& beta, const number& a, const number& f, mat& m) {
	number res = 0.0;
	
	uint mj = negNeighDisjoint(j,l.size());
	uint pj = posNeighDisjoint(j,l.size());	
	uint mk = negNeighDisjoint(k,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	
	number x_jk = DistanceDisjoint(l[j],l[k],beta);
	number x_jmk = DistanceDisjoint(l[j],l[mk],beta);
	number x_mjk = DistanceDisjoint(l[mj],l[k],beta);
	number x_mjmk = DistanceDisjoint(l[mj],l[mk],beta);
	
	number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0*FVacuum(x_mjmk,a) \
			 + 2.0*FVacuum(x_jk,a) \
			 - 2.0*dFVacuumdronr(x_jk,a)*T_jk;
		if (k!=mj)
			res += - 2.0*FVacuum(x_mjk,a);
		if (k!=pj)
			res += - 2.0*FVacuum(x_jmk,a);
	}

	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res +=  + 2.0*dFVacuumdronr(x_jk,a)*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,k,mu,beta) \
				 - 2.0*dFVacuumdronr(x_jk,a)*DXDisjoint(l,j,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta) \
				 - (2.0*ddFVacuumdrdr(x_jk,a)*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(x_jk,2) \
				 + (2.0*dFVacuumdronr(x_jk,a)*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(x_jk,2);
	if (k!=mj)
		res += + 2.0*dFVacuumdronr(x_mjk,a)*DXDisjoint(l,mj,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta);
	if (k!=pj)
		res += - 2.0*dFVacuumdronr(x_jmk,a)*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,mk,mu,beta);
		
	// terms with sums
	if (k==j || k==mj || k==pj) {
	
		uint pi;
		number x_ij, x_imj, T_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			pi = (i==(l.size()-1)? 0: i+1);
			x_ij = DistanceDisjoint(l[i],l[j],beta);
			x_imj = DistanceDisjoint(l[i],l[mj],beta);
			T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
			
			if (k==j && i!=j) {
				res +=  2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,j,i,nu,beta) \
						 + 2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta) \
						 + (2.0*ddFVacuumdrdr(x_ij,a)*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(x_ij,2) \
						 - (2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(x_ij,2);
				if (mu==nu)
					res += + 2.0*dFVacuumdronr(x_ij,a)*T_ij; //
			}
			if (k==mj && i!=mj) 
				res += 2.0*dFVacuumdronr(x_imj,a)*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,i,k,nu,beta);
			if (k==pj && i!=j) 
				res += - 2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta);
			
		}		
	}
	
	//coincident terms
	if (k==j && mu==nu)
		res += (2.0*FVacuum(0.0,a));
	if (k==mj && mu==nu)
		res += (-1.0*FVacuum(0.0,a));
	if (k==pj && mu==nu)
		res += (-1.0*FVacuum(0.0,a));
	
	m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
	
}

// mdV0DisjointLRGeneric
template<uint Dim>
static void mdV0DisjointLRGeneric(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					VacuumFunction FVacuum, VacuumFunction dFVacuumdronr, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && i<N/2) || (j>=N/2 && i>=N/2))
		mdV0DisjointGeneric(j,mu,i,l,FVacuum,dFVacuumdronr,beta,a,f,v);
}

// ddV0DisjointLRGeneric
template<uint Dim>
static void ddV0DisjointLRGeneric(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						VacuumFunction FVacuum, VacuumFunction dFVacuumdronr,VacuumFunction ddFVacuumdrdr,\
						 const number& beta, const number& a, const number& f, mat& m) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2)) {
		number res = 0.0;
	
		uint mj = negNeighDisjoint(j,l.size());
		uint pj = posNeighDisjoint(j,l.size());	
		uint mk = negNeighDisjoint(k,l.size());
		uint pk = posNeighDisjoint(k,l.size());
	
		number x_jk = DistanceDisjoint(l[j],l[k],beta);
		number x_jmk = DistanceDisjoint(l[j],l[mk],beta);
		number x_mjk = DistanceDisjoint(l[mj],l[k],beta);
		number x_mjmk = DistanceDisjoint(l[mj],l[mk],beta);
	
		number T_jk = DotDisjoint(l[pj],l[j],l[pk],l[k],beta);
	
		// terms where mu==nu, without sums
		if (mu==nu) {
			if (k!=j)
				res += 2.0*FVacuum(x_mjmk,a) \
				 + 2.0*FVacuum(x_jk,a) \
				 - 2.0*dFVacuumdronr(x_jk,a)*T_jk;
			if (k!=mj)
				res += - 2.0*FVacuum(x_mjk,a);
			if (k!=pj)
				res += - 2.0*FVacuum(x_jmk,a);
		}

		// terms where mu not nexcessarily equal to nu, without sums
		if (k!=j)
			res +=  + 2.0*dFVacuumdronr(x_jk,a)*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,k,mu,beta) \
					 - 2.0*dFVacuumdronr(x_jk,a)*DXDisjoint(l,j,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta) \
					 - (2.0*ddFVacuumdrdr(x_jk,a)*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(x_jk,2) \
					 + (2.0*dFVacuumdronr(x_jk,a)*DXDisjoint(l,j,k,mu,beta)*DXDisjoint(l,j,k,nu,beta)*T_jk)/pow(x_jk,2);
		if (k!=mj)
			res += + 2.0*dFVacuumdronr(x_mjk,a)*DXDisjoint(l,mj,k,nu,beta)*DXDisjoint(l,k,pk,mu,beta);
		if (k!=pj)
			res += - 2.0*dFVacuumdronr(x_jmk,a)*DXDisjoint(l,j,pj,nu,beta)*DXDisjoint(l,j,mk,mu,beta);
		
		// terms with sums
		if (k==j || k==mj || k==pj) {
	
			uint pi;
			number x_ij, x_imj, T_ij;
		
			// for LR
			uint iMin = (j<N/2? 0: N/2);
			uint iMax = (j<N/2? N/2:N);
		
			for (uint i=iMin; i<iMax; i++) {
		
				pi = (i==(l.size()-1)? 0: i+1);
				x_ij = DistanceDisjoint(l[i],l[j],beta);
				x_imj = DistanceDisjoint(l[i],l[mj],beta);
				T_ij = DotDisjoint(l[pi],l[i],l[pj],l[j],beta);
			
				if (k==j && i!=j) {
					res +=  2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,j,i,nu,beta) \
							 + 2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta) \
							 + (2.0*ddFVacuumdrdr(x_ij,a)*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(x_ij,2) \
							 - (2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,j,i,mu,beta)*DXDisjoint(l,j,i,nu,beta)*T_ij)/pow(x_ij,2);
					if (mu==nu)
						res += + 2.0*dFVacuumdronr(x_ij,a)*T_ij; //
				}
				if (k==mj && i!=mj) 
					res += 2.0*dFVacuumdronr(x_imj,a)*DXDisjoint(l,i,pi,mu,beta)*DXDisjoint(l,i,k,nu,beta);
				if (k==pj && i!=j) 
					res += - 2.0*dFVacuumdronr(x_ij,a)*DXDisjoint(l,i,pi,nu,beta)*DXDisjoint(l,j,i,mu,beta);
			
			}		
		}
	
		//coincident terms
		if (k==j && mu==nu)
			res += (2.0*FVacuum(0.0,a));
		if (k==mj && mu==nu)
			res += (-1.0*FVacuum(0.0,a));
		if (k==pj && mu==nu)
			res += (-1.0*FVacuum(0.0,a));
	
		m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
	}	
}


/*----------------------------------------------------------------------------------------------------------------------------
	1.0 - nr loop scalar functions
----------------------------------------------------------------------------------------------------------------------------*/

// L
template <uint Dim>
void NEWL (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pos = posNeigh(j,l.size());
	result += f*Distance(l[pos],l[j]);
}

// LDisjoint
template <uint Dim>
void NEWLDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pos = posNeighDisjoint(j,l.size());
	result += f*DistanceDisjoint(l[pos],l[j],beta);
}

// DistPow
template <uint Dim>
void NEWDistPow (const uint& j, const Loop<Dim>& l, const number& w, const number& f, number& result) {
	uint pos = posNeigh(j,l.size());
	result += f*pow(DistanceSquared(l[pos],l[j]),w/2.0);
}

// Sm
template <uint Dim>
void NEWSm (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
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
void NEWS0 (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	result += f*DistanceSquared(l[pj],l[j])*(number)l.size()/4.0;
}


// S0Disjoint
template <uint Dim>
void NEWS0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	result += f*DistanceSquaredDisjoint(l[pj],l[j],beta)*(number)l.size()/4.0;
}

// I0
template <uint Dim>
void NEWI0 (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	result += f*(l[j])[Dim-2]*DX(l,pj,j,Dim-1);
}

// I0Disjoint
template <uint Dim>
void NEWI0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	result += f*(l[j])[Dim-2]*DXDisjoint(l,pj,j,Dim-1,beta);
}

// In
template <uint Dim>
void NEWIn (const uint& j, const Loop<Dim>& l, const number& n, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	result += f*pow((l[j])[Dim-2],n)*DX(l,pj,j,Dim-1);
}

// InDisjoint
template <uint Dim>
void NEWInDisjoint (const uint& j, const Loop<Dim>& l, const number& n, const number& beta, const number& f, number& result) {
	uint pj = posNeighDisjoint(j,l.size());
	result += f*pow((l[j])[Dim-2],n)*DXDisjoint(l,pj,j,Dim-1,beta);
}

// Angle
template <uint Dim>
void NEWAngle (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	number cos_gamma;
	
	cos_gamma = Dot(l[pj],l[j],l[j],l[nj]);
	cos_gamma /= Distance(l[pj],l[j])*Distance(l[j],l[nj]);
	
	result += f*(acos(cos_gamma));
}

// NEWPseudoAngle, in 3-4 plane
template <uint Dim>
void NEWPseudoAngle (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	number sin_gamma;
	
	sin_gamma = (l[pj]-l[j])[Dim-2]*(l[j]-l[nj])[Dim-1] - (l[pj]-l[j])[Dim-1]*(l[j]-l[nj])[Dim-2];
	sin_gamma /= sqrt((pow((l[pj]-l[j])[Dim-2],2)+pow((l[pj]-l[j])[Dim-1],2))*(pow((l[j]-l[nj])[Dim-2],2)+pow((l[j]-l[nj])[Dim-1],2)));
	
	result += f*(asin(sin_gamma));
}

// FGamma
template <uint Dim>
void NEWFGamma (const uint& j,const Loop<Dim>& l,  const number& f, number& result) {

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
void NEWMaxXn (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result) {
	if (n>=Dim) {
		cerr << "NEWMaxXn error: " << n << ">=" << Dim << endl;
	}
	if (k<j) {
		number temp = abs((l[j])[n]-(l[k])[n]);
		result = (temp*f>result? temp*f: result);
	}
}

// MinXnDisjoint
template <uint Dim>
void NEWMinXnDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result) {
	if (n>=Dim) {
		cerr << "NEWMinXnDisjoint error: " << n << ">=" << Dim << endl;
	}
	if (k<(l.size()/2) && j>=(l.size()/2)) {
		number temp = abs((l[j])[n]-(l[k])[n]);
		result = (temp*f<result? temp*f: result);
	}
}

// Vor
template <uint Dim>
void NEWVor (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {
	V0Generic(j,k,l,&F0,a,f,result);
}

// Vlr
template <uint Dim>
void NEWVlr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {
	if (k<=j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*pow(0.25*DistanceSquared(l[pj]+l[j],l[pk]+l[k])+a*a,(2.0-Dim)/2.0);	
	}
}

// Ver
template <uint Dim>
void NEWVer (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {
	V0Generic(j,k,l,&FExp,a,f,result);
}

// Vdr
template <uint Dim>
void NEWVdr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {
	V0Generic(j,k,l,&FDim,a,f,result);
}

// Vthr
template <uint Dim>
void NEWVthr (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	VthGeneric(j,k,l,&FThermal,beta,a,f,result);
}

// VthrDisjoint
template <uint Dim>
void NEWVthrDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	VthDisjointGeneric(j,k,l,&FThermal,beta,a,f,result);
}

// VthrDisjointLR
template <uint Dim>
void NEWVthrDisjointLR (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	VthDisjointLR2Generic(j,k,l,&FThermal,beta,a,f,result);
}

// VnonrelDisjoint
template <uint Dim>
void NEWVnonrelDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
	if (j<l.size()/2) {
		number dt = 2.0*beta/(number)l.size();
		uint oj = oppNeigh(j,l.size());

		result += f*dt/SpatialDistance(l[j],l[oj]); // dt/r, not lorentz invariant but galilean invariant
	}
}

// VnonrelrDisjoint
template <uint Dim>
void NEWVnonrelrDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	if (j<l.size()/2) {
		uint oj = oppNeigh(j,l.size());
		number dt = 2.0*beta/(number)l.size();
		number r2 = SpatialDistanceSquared(l[j],l[oj]);
		number a2 = pow(a,2);

		result += f*dt*(1.0/sqrt(r2 + a2) - exp(-r2/a2)/a); // not lorentz invariant but galilean invariant
	}
}

// Gaussian
template <uint Dim>
void NEWGaussian (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {
	V0Generic(j,k,l,&G0,a,f,result);
}

// GaussianDisjoint
template <uint Dim>
void NEWGaussianDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	V0DisjointGeneric(j,k,l,&G0,beta,a,f,result);
}

// GaussianLRDisjoint
template <uint Dim>
void NEWGaussianLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	V0DisjointLRGeneric(j,k,l,&G0,beta,a,f,result);
}

// GaussianThermal
template <uint Dim>
void NEWGaussianThermal (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	VthGeneric(j,k,l,&GThermal,beta,a,f,result);
}

// NEWGaussianThermal2
template <uint Dim>
void NEWGaussianThermal2 (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result) {
	VthGeneric(j,k,l,&G2Thermal,beta,a,f,result);
}

// GaussianThermalDisjoint
template <uint Dim>
void NEWGaussianThermalDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	VthDisjointGeneric(j,k,l,&GThermal,beta,a,f,result);
}

// NEWGaussianThermal2Disjoint
template <uint Dim>
void NEWGaussianThermal2Disjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	VthDisjointGeneric(j,k,l,&G2Thermal,beta,a,f,result);
}

// GaussianThermalLRDisjoint
template <uint Dim>
void NEWGaussianThermalLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	VthDisjointLRGeneric(j,k,l,&GThermal,beta,a,f,result);
}

// NEWGaussianThermal2LRDisjoint
template <uint Dim>
void NEWGaussianThermal2LRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result) {
	VthDisjointLRGeneric(j,k,l,&G2Thermal,beta,a,f,result);
}

// InlineCurvatureMax
template <uint Dim>
void NEWInlineCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	number temp = f*abs( Dot(l[pj],l[j],l[pj],l[j]) + Dot(l[nj],l[j],l[pj],l[j]) )/DistanceSquared(l[pj],l[j]);
	result = (temp>result? temp: result);
}

// InlineCurvatureMax
template <uint Dim>
void NEWInlineCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
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
void NEWInlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result){
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	number temp = f*abs( DotDisjoint(l[pj],l[j],l[pj],l[j],beta) + DotDisjoint(l[nj],l[j],l[pj],l[j],beta) )/DistanceSquaredDisjoint(l[pj],l[j],beta);
	result = (temp>result? temp: result);
}

// InlineCurvatureMaxDisjoint
template <uint Dim>
void NEWInlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
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
void NEWCuspCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
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
void NEWCuspCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
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
void NEWCuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
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
void NEWCuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
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
void NEWKGMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
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
void NEWKGMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
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
void NEWKGMaxPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
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
void NEWKGMaxPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
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
void NEWKGMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
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
void NEWKGMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
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
void NEWKGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result) {
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
void NEWKGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result) {
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


/*----------------------------------------------------------------------------------------------------------------------------
	1.1 - nr loop vector and matrix functions
----------------------------------------------------------------------------------------------------------------------------*/

// mdPX_nr
template<uint Dim>
void NEWmdPX_nr(const Loop<Dim>& l,const uint& loc, const Point<Dim>& P, const number& f, vec& v) {
	for (uint mu=0; mu<Dim; mu++)
		v[loc*Dim+mu] += -f*P[mu];
}

// mdL_nr
template<uint Dim>
void NEWmdL_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
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
void NEWmdLDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
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
void NEWddL_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
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
void NEWddLDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
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
void NEWmdDistPow_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& w,  const number& f, vec& v) {
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
void NEWddDistPow_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& w\
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
void NEWmdI0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
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
void NEWmdI0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
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
void NEWddI0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
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
void NEWddI0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& beta, const number& f, mat& m) {
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
void NEWmdIn_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& f, vec& v) {
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
void NEWmdInDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, vec& v) {
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
void NEWddIn_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& f, mat& m) {
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
void NEWddInDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, mat& m) {
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
void NEWmdS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	v[j*Dim+mu] += -f*(DX(l,j,nj,mu) + DX(l,j,pj,mu))*(number)l.size()/2.0;
}

// mdS0Disjoint_nr
template<uint Dim>
void NEWmdS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	v[j*Dim+mu] += -f*(DXDisjoint(l,j,nj,mu,beta) + DXDisjoint(l,j,pj,mu,beta))*(number)l.size()/2.0;
}

// ddS0_nr
template<uint Dim>
void NEWddS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
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
void NEWddS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
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
void NEWmdsqrtS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0, const number& f, vec& v) {
	NEWmdS0_nr(j,mu,l,f*2.0/sqrt4s0,v);
}

// mdsqrtS0Disjoint_nr
template<uint Dim>
void NEWmdsqrtS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0,\
 const number& beta, const number& f, vec& v) {
 	NEWmdS0Disjoint_nr(j,mu,l,beta,f*2.0/sqrt4s0,v);
 }

// ddsqrtS0_nr
template<uint Dim>
void NEWddsqrtS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& f, mat& m) {
	NEWddS0_nr(j,mu,k,nu,l,f*2.0/sqrt4s0,m);
	
	uint pj = posNeigh(j,l.size());
	uint nj = negNeigh(j,l.size());
	uint pk = posNeigh(k,l.size());
	uint nk = negNeigh(k,l.size());
	
	m(j*Dim+mu,k*Dim+nu) -= (f*pow((number)l.size(),2)/pow(sqrt4s0,3)) * (DX(l,j,nj,mu) + DX(l,j,pj,mu)) \
												* (DX(l,k,nk,nu) + DX(l,k,pk,nu));
}


// ddsqrtS0Disjoint_nr
template<uint Dim>
void NEWddsqrtS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& beta, const number& f, mat& m) {
	NEWddS0Disjoint_nr(j,mu,k,nu,l,f*2.0/sqrt4s0,m);
	
	uint pj = posNeighDisjoint(j,l.size());
	uint nj = negNeighDisjoint(j,l.size());
	uint pk = posNeighDisjoint(k,l.size());
	uint nk = negNeighDisjoint(k,l.size());
	
	m(j*Dim+mu,k*Dim+nu) -= (f*pow((number)l.size(),2)/pow(sqrt4s0,3)) * (DXDisjoint(l,j,nj,mu,beta) + DXDisjoint(l,j,pj,mu,beta)) \
												* (DXDisjoint(l,k,nk,nu,beta) + DXDisjoint(l,k,pk,nu,beta));				
}

// mdVor_nr
template<uint Dim>
void NEWmdVor_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	mdV0Generic(j, mu, i, l, &F0, &DF0DrOnr, a, f, v);
}

// ddVor_nr
template<uint Dim>
void NEWddVor_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	ddV0Generic(j, mu, k, nu, l, &F0, &DF0DrOnr, &DDF0DrDr, a, f, m);
}

// mdVlr_nr
template<uint Dim>
void NEWmdVlr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
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

	v[j*Dim+mu] += -f*res;
}

// ddVlr_nr
template<uint Dim>
void NEWddVlr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
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
	
	m(Dim*j+mu,Dim*k+nu) += f*res;
}

// mdVer_nr
template<uint Dim>
void NEWmdVer_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	mdV0Generic(j, mu, i, l, &FExp, &DFExpDrOnr, a, f, v);
}

// ddVer_nr
template<uint Dim>
void NEWddVer_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	ddV0Generic(j, mu, k, nu, l, &FExp, &DFExpDrOnr, &DDFExpDrDr, a, f, m);
}

// mdVdr_nr
template<uint Dim>
void NEWmdVdr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	mdV0Generic(j, mu, i, l, &FDim, &DFDimDrOnr, a, f, v);
}

// ddVdr_nr
template<uint Dim>
void NEWddVdr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	ddV0Generic(j, mu, k, nu, l, &FDim, &DFDimDrOnr, &DDFDimDrDr, a, f, m);
}

// mdVnonrelDisjoint_nr
template<uint Dim>
void NEWmdVnonrelDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v) {
	if (mu<(Dim-1)) {
		uint oj = oppNeigh(j,l.size());
		number dt = 2.0*beta/(number)l.size();
		v[j*Dim+mu] += -f*dt*DX(l,j,oj,mu)*(-1.0/pow(SpatialDistance(l[j],l[oj]),3));
	}
}

// mdVnonrelrDisjoint_nr
template<uint Dim>
void NEWmdVnonrelrDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& a\
						, const number& f, vec& v) {
	if (mu<(Dim-1)) {
		uint oj = oppNeigh(j,l.size());
		number dt = 2.0*beta/(number)l.size();
		number r2 = SpatialDistanceSquared(l[j],l[oj]);
		number a2 = pow(a,2);
		v[j*Dim+mu] += -f*dt*DX(l,j,oj,mu)*(-1.0/pow(r2 + a2,3.0/2.0) + 2.0*exp(-r2/a2)/a/a2 );
	}
}

// mdVthr_nr
template<uint Dim>
void NEWmdVthr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v) {
	mdVthGeneric(j, mu, i, l, &FThermal, &DFThermalDrOnr,&DFThermalDt, beta, a, f, v);
}

// ddVthr_nr
template<uint Dim>
void NEWddVthr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthGeneric(j, mu, k, nu, l, &FThermal, &DFThermalDrOnr,&DFThermalDt,&DDFThermalDrDr,&DDFThermalDtDrOnr,&DDFThermalDtDt, beta, a, f, m);
}

// mdVthrDisjoint_nr
template<uint Dim>
void NEWmdVthrDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v) {
	mdVthDisjointGeneric(j, mu, i, l, &FThermal, &DFThermalDrOnr,&DFThermalDt, beta, a, f, v);
}

// mdVthrDisjointLR_nr
template<uint Dim>
void NEWmdVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v) {
	mdVthDisjointLR2Generic(j, mu, i, l, &FThermal, &DFThermalDrOnr,&DFThermalDt, beta, a, f, v);
}

// ddVthrDisjoint_nr
template<uint Dim>
void NEWddVthrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthDisjointGeneric(j, mu, k, nu, l, &FThermal, &DFThermalDrOnr,&DFThermalDt,&DDFThermalDrDr,&DDFThermalDtDrOnr,&DDFThermalDtDt, beta, a, f, m);
}

// ddVthrDisjointLR_nr
template<uint Dim>
void NEWddVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthDisjointLR2Generic(j, mu, k, nu, l, &FThermal, &DFThermalDrOnr,&DFThermalDt,&DDFThermalDrDr,&DDFThermalDtDrOnr,&DDFThermalDtDt, beta, a, f, m);
}


// ddVnonrelDisjoint_nr
template<uint Dim>
void NEWddVnonrelDisjoint_nr (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& f, mat& m) {
	
	if (mu!=(Dim-1) && nu!=(Dim-1)) {
		uint oj = oppNeigh(j,l.size());
		uint ok = oppNeigh(k,l.size());
		number dt = 2.0*beta/(number)l.size();
	
		if (j==k) {
			m(Dim*j+mu,Dim*k+nu) += -f*dt*(-3.0)*DX(l,j,oj,mu)*DX(l,j,oj,nu)\
								/pow(SpatialDistance(l[j],l[oj]),5);
			if (mu==nu)
				m(Dim*j+mu,Dim*k+nu) += -f*dt/pow(SpatialDistance(l[j],l[oj]),3);
		}
		else if (k==oj) {
			m(Dim*j+mu,Dim*k+nu) += -f*dt*(-3.0)*DX(l,j,oj,mu)*DX(l,k,ok,nu)\
								/pow(SpatialDistance(l[j],l[oj]),5);
			if (mu==nu)
				m(Dim*j+mu,Dim*k+nu) += f*dt/pow(SpatialDistance(l[j],l[oj]),3);
		}
	}
	
}

// ddVnonrelrDisjoint_nr
template<uint Dim>
void NEWddVnonrelrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {

	if (mu!=(Dim-1) && nu!=(Dim-1)) {
		uint oj = oppNeigh(j,l.size());
		uint ok = oppNeigh(k,l.size());
		number dt = 2.0*beta/(number)l.size();
		number r2 = SpatialDistanceSquared(l[j],l[oj]);
		number a2 = pow(a,2);	
	
		if (j==k) {
			m(Dim*j+mu,Dim*k+nu) += f*dt*DX(l,j,oj,mu)*DX(l,j,oj,nu)\
								*( 3.0/pow(r2+a2,5.0/2.0) - 4.0*exp(-r2/a2)/pow(a,5) );
			if (mu==nu)
				m(Dim*j+mu,Dim*k+nu) += f*dt*(-1.0/pow(r2 + a2,3.0/2.0) + 2.0*exp(-r2/a2)/pow(a,3) );
		}
		else if (k==oj) {
			m(Dim*j+mu,Dim*k+nu) += f*dt*DX(l,j,oj,mu)*DX(l,k,ok,nu)\
								*( 3.0/pow(r2+a2,5.0/2.0) - 4.0*exp(-r2/a2)/pow(a,5) );
			if (mu==nu)
				m(Dim*j+mu,Dim*k+nu) += -f*dt*(-1.0/pow(r2 + a2,3.0/2.0) + 2.0*exp(-r2/a2)/pow(a,3) );
		}
	}
	
}

// mdGaussian_nr
template<uint Dim>
void NEWmdGaussian_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	mdV0Generic(j, mu, i, l, &G0, &DG0DrOnr, a, f, v);
}

// mdGaussianThermal_nr
template<uint Dim>
void NEWmdGaussianThermal_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v) {
	mdVthGeneric(j, mu, i, l, &GThermal, &DGThermalDrOnr,&DGThermalDt, beta, a, f, v);
}

// mdGaussianThermal2_nr
template<uint Dim>
void NEWmdGaussianThermal2_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v) {
	mdVthGeneric(j, mu, i, l, &G2Thermal, &DG2ThermalDrOnr,&DG2ThermalDt, beta, a, f, v);
}

// ddGaussian_nr
template<uint Dim>
void NEWddGaussian_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	ddV0Generic(j, mu, k, nu, l, &G0, &DG0DrOnr, &DDG0DrDr, a, f, m);
}


// ddGaussianThermal_nr
template<uint Dim>
void NEWddGaussianThermal_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthGeneric(j, mu, k, nu, l, &GThermal, &DGThermalDrOnr,&DGThermalDt,&DDGThermalDrDr,&DDGThermalDtDrOnr,&DDGThermalDtDt, beta, a, f, m);
}

// ddGaussianThermal2_nr
template<uint Dim>
void NEWddGaussianThermal2_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthGeneric(j, mu, k, nu, l, &G2Thermal, &DG2ThermalDrOnr,&DG2ThermalDt,&DDG2ThermalDrDr,&DDG2ThermalDtDrOnr,&DDG2ThermalDtDt, beta, a, f, m);
}

// mdGaussianDisjoint_nr
template<uint Dim>
void NEWmdGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	mdV0DisjointGeneric(j, mu, i, l, &G0, &DG0DrOnr, beta, a, f, v);
}

// mdGaussianLRDisjoint_nr
template<uint Dim>
void NEWmdGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	mdV0DisjointLRGeneric(j, mu, i, l, &G0, &DG0DrOnr, beta, a, f, v);
}


// mdGaussianThermalDisjoint_nr
template<uint Dim>
void NEWmdGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	mdVthDisjointGeneric(j, mu, i, l, &GThermal, &DGThermalDrOnr,&DGThermalDt, beta, a, f, v);
}

// mdGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWmdGaussianThermal2Disjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	mdVthDisjointGeneric(j, mu, i, l, &G2Thermal, &DG2ThermalDrOnr,&DG2ThermalDt, beta, a, f, v);
}
			
// mdGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWmdGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	mdVthDisjointLRGeneric(j, mu, i, l, &GThermal, &DGThermalDrOnr,&DGThermalDt, beta, a, f, v);
}

// mdGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWmdGaussianThermal2LRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v) {
	mdVthDisjointLRGeneric(j, mu, i, l, &G2Thermal, &DG2ThermalDrOnr,&DG2ThermalDt, beta, a, f, v);
}
		 
// ddGaussianDisjoint_nr
template<uint Dim>
void NEWddGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddV0DisjointGeneric(j, mu, k, nu, l, &G0, &DG0DrOnr, &DDG0DrDr, beta, a, f, m);
}

// ddGaussianLRDisjoint_nr
template<uint Dim>
void NEWddGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddV0DisjointLRGeneric(j, mu, k, nu, l, &G0, &DG0DrOnr, &DDG0DrDr, beta, a, f, m);
}
						 
// ddGaussianThermalDisjoint_nr
template<uint Dim>
void NEWddGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthDisjointGeneric(j, mu, k, nu, l, &GThermal, &DGThermalDrOnr,&DGThermalDt,&DDGThermalDrDr,&DDGThermalDtDrOnr,&DDGThermalDtDt, beta, a, f, m);
}
 
// ddGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWddGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthDisjointLRGeneric(j, mu, k, nu, l, &GThermal, &DGThermalDrOnr,&DGThermalDt,&DDGThermalDrDr,&DDGThermalDtDrOnr,&DDGThermalDtDt, beta, a, f, m);	 
}

// ddGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWddGaussianThermal2Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthDisjointGeneric(j, mu, k, nu, l, &G2Thermal, &DG2ThermalDrOnr,&DG2ThermalDt,&DDG2ThermalDrDr,&DDG2ThermalDtDrOnr,&DDG2ThermalDtDt, beta, a, f, m);
}
 
// ddGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWddGaussianThermal2LRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m) {
	ddVthDisjointLRGeneric(j, mu, k, nu, l, &G2Thermal, &DG2ThermalDrOnr,&DG2ThermalDt,&DDG2ThermalDrDr,&DDG2ThermalDtDrOnr,&DDG2ThermalDtDt, beta, a, f, m);	 
}

// mdGamma_nr
template<uint Dim>
void NEWmdFGamma_nr(const Loop<Dim>& l, const uint& j, const number& f, vec& v) {
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
void NEWddFGamma_nr(const Loop<Dim>& l, const uint& j, const number& f, mat& m) {
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
void NEWErgS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, number& erg) {
	erg += f*DX(l,j,mu)*(number)l.size()/2.0;
}

// ErgS0Disjoint_nr
template<uint Dim>
void NEWErgS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, number& erg) {
	erg += f*DXDisjoint(l,j,mu,beta)*(number)l.size()/2.0;
}

// ErgL_nr
template<uint Dim>
void NEWErgL_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, number& erg) {
	uint pj = posNeigh(j,l.size());
	erg += f*DX(l,pj,j,mu)/Distance(l[pj],l[j]);
}

// ErgLDisjoint_nr
template<uint Dim>
void NEWErgLDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, number& erg) {
	uint pj = posNeighDisjoint(j,l.size());
	erg += f*DXDisjoint(l,pj,j,mu,beta)/DistanceDisjoint(l[pj],l[j],beta);
}

// ErgsqrtS0_nr
template<uint Dim>
void NEWErgsqrtS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& f, number& erg) {
	NEWErgS0_nr(l,j,mu,f*2.0/sqrt4s0,erg);
}

// ErgsqrtS0Disjoint_nr
template<uint Dim>
void NEWErgsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, number& erg) {
	NEWErgS0Disjoint_nr(l,j,mu,beta,f*2.0/sqrt4s0,erg);
}

// ErgI0_nr
template<uint Dim>
void NEWErgI0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeigh(j,l.size());
		erg += -f*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*(l[j])[Dim-2];
}

// ErgI0Disjoint_nr
template<uint Dim>
void NEWErgI0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeighDisjoint(j,l.size());
		erg += -f*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*(l[j])[Dim-2];
}

// ErgIn_nr
template<uint Dim>
void NEWErgIn_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeigh(j,l.size());
		erg += -f*n*pow((l[j])[Dim-2],n-1.0)*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*pow((l[j])[Dim-2],n);
}

// ErgInDisjoint_nr
template<uint Dim>
void NEWErgInDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& beta, const number& f, number& erg) {
	if (mu==(Dim-2)) {
		uint pj = posNeighDisjoint(j,l.size());
		erg += -f*n*pow((l[j])[Dim-2],n-1.0)*(l[pj])[Dim-1];
	}
	else if (mu==(Dim-1))
		erg += f*pow((l[j])[Dim-2],n);
}

// ErgVor_nr
template<uint Dim>
void NEWErgVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg) {
	erg += f*DX(l,k,mu)/(DistanceSquared(l[j],l[k]) + a*a);
}

// ErgVthr_nr
template<uint Dim>
void NEWErgVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);

	erg += f*(-pow(2.0*PI,2))*DX(l,k,mu)*FThermal(r,t,beta,a);
}

// ErgVthrDisjoint_nr
template<uint Dim>
void NEWErgVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1); // intentionally not DXDisjoint(l,k,j,Dim-1,beta);

	erg += f*(-pow(2.0*PI,2))*DXDisjoint(l,k,mu,beta)*FThermal(r,t,beta,a);	
}

// ErgVthrDisjointLR_nr
template<uint Dim>
void NEWErgVthrDisjointLR_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWErgVthrDisjoint_nr(l,j,mu,k,beta,a,f,erg);
	else
		NEWErgVthrDisjoint_nr(l,j,mu,k,beta,0.0,f,erg);
}

// ErgGaussian_nr
template<uint Dim>
void NEWErgGaussian_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg) {
	erg += f*DX(l,k,mu)*exp(-DistanceSquared(l[j],l[k])/a/a);
}

// ErgGaussianDisjoint_nr
template<uint Dim>
void NEWErgGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	erg += f*DXDisjoint(l,k,mu,beta)*exp(-DistanceSquared(l[j],l[k])/a/a); //intentionally not disjoint
}

// ErgGaussianLRDisjoint_nr
template<uint Dim>
void NEWErgGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWErgGaussianDisjoint_nr(l,j,mu,k,beta,a,f,erg);
}

// ErgGaussianThermal_nr
template<uint Dim>
void NEWErgGaussianThermal_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DX(l,k,mu)*GThermal(r,t,beta,a);
}

// ErgGaussianThermal2_nr
template<uint Dim>
void NEWErgGaussianThermal2_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DX(l,k,mu)*G2Thermal(r,t,beta,a);
}

// ErgGaussianThermalDisjoint_nr
template<uint Dim>
void NEWErgGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DXDisjoint(l,k,mu,beta)*GThermal(r,t,beta,a);
}

// ErgGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWErgGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	number r = SpatialDistance(l[j],l[k]);
	number t = DX(l,k,j,Dim-1);
	erg += f*DXDisjoint(l,k,mu,beta)*G2Thermal(r,t,beta,a);
}

// ErgGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWErgGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWErgGaussianThermalDisjoint_nr(l,j,mu,k,beta,a,f,erg);
}

// ErgGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWErgGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWErgGaussianThermal2Disjoint_nr(l,j,mu,k,beta,a,f,erg);
}

// PS0_nr
template<uint Dim>
void NEWPS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, vec& v) {
	NEWErgS0_nr(l,j,mu,f,v[Dim*j+mu]);
}

// PS0Disjoint_nr
template<uint Dim>
void NEWPS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, vec& v) {
	NEWErgS0Disjoint_nr(l,j,mu,beta,f,v[Dim*j+mu]);
}

// PL_nr
template<uint Dim>
void NEWPL_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, vec& v) {
	NEWErgL_nr(l,j,mu,f,v[Dim*j+mu]);
}

// PLDisjoint_nr
template<uint Dim>
void NEWPLDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, vec& v) {
	NEWErgLDisjoint_nr(l,j,mu,beta,f,v[Dim*j+mu]);
}

// PsqrtS0_nr
template<uint Dim>
void NEWPsqrtS0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& f, vec& v) {
	NEWPS0_nr(l,j,mu,f*2.0/sqrt4s0,v);
}

// PsqrtS0Disjoint_nr
template<uint Dim>
void NEWPsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, vec& v) {
	NEWPS0Disjoint_nr(l,j,mu,beta,f*2.0/sqrt4s0,v);
}

// PI0_nr
template<uint Dim>
void NEWPI0_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& f, vec& v) {
	NEWErgI0_nr(l,j,mu,f,v[Dim*j+mu]);
}

// PI0Disjoint_nr
template<uint Dim>
void NEWPI0Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& beta, const number& f, vec& v) {
	NEWErgI0Disjoint_nr(l,j,mu,beta,f,v[Dim*j+mu]);
}

// PIn_nr
template<uint Dim>
void NEWPIn_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& f, vec& v) {
	NEWErgIn_nr(l,j,mu,n,f,v[Dim*j+mu]);
}

// PInDisjoint_nr
template<uint Dim>
void NEWPInDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const number& n, const number& beta, const number& f, vec& v) {
	NEWErgInDisjoint_nr(l,j,mu,n,beta,f,v[Dim*j+mu]);
}

// PVor_nr
template<uint Dim>
void NEWPVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v) {
	NEWErgVor_nr(l,j,mu,k,a,2.0*f,v[Dim*j+mu]);
}

// PRVor_nr
template<uint Dim>
void NEWPRVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v) {
	number res = 0.0;
	for (uint i=0; i<=j; i++) {
		res += DX(l,k,nu)*DX(l,i,nu)*(-2.0*DX(l,k,i,mu))/pow(DistanceSquared(l[i],l[k]) + a*a,2);		
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PVthr_nr
template<uint Dim>
void NEWPVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {	
	NEWErgVthr_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRVthr_nr
template<uint Dim>
void NEWPRVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
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
void NEWPVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {	
	NEWErgVthrDisjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRVthrDisjoint_nr
template<uint Dim>
void NEWPRVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
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
void NEWPGaussian_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v) {
	NEWErgGaussian_nr(l,j,mu,k,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussian_nr
template<uint Dim>
void NEWPRGaussian_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v) {
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
void NEWPGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	NEWErgGaussianDisjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussianDisjoint_nr
template<uint Dim>
void NEWPRGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
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
void NEWPGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWPGaussianDisjoint_nr(l,j,mu,k,beta,a,f,v);
}

// PRGaussianLRDisjoint_nr
template<uint Dim>
void NEWPRGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWPRGaussianDisjoint_nr(l,j,mu,k,nu,beta,a,f,v);
}

// PGaussianThermal_nr
template<uint Dim>
void NEWPGaussianThermal_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	NEWErgGaussianThermal_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PGaussianThermal2_nr
template<uint Dim>
void NEWPGaussianThermal2_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	NEWErgGaussianThermal2_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussianThermal_nr
template<uint Dim>
void NEWPRGaussianThermal_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
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
void NEWPRGaussianThermal2_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeigh(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += DX(l,pk,k,nu)*DX(l,pi,i,nu)*DX(l,k,i,mu)*DG2ThermalDrOnr(r,t,beta,a);
		}
		else {
			res += DX(l,pk,k,nu)*DX(l,pi,i,nu)*DG2ThermalDt(r,t,beta,a);
		}	
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PGaussianThermalDisjoint_nr
template<uint Dim>
void NEWPGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	NEWErgGaussianThermalDisjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWPGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	NEWErgGaussianThermal2Disjoint_nr(l,j,mu,k,beta,a,2.0*f,v[Dim*j+mu]);
}

// PRGaussianThermalDisjoint_nr
template<uint Dim>
void NEWPRGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
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
void NEWPRGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint pk = posNeighDisjoint(k,l.size()), pi;
	number res = 0.0, r, t;
	for (uint i=0; i<=j; i++) {
		pi = posNeighDisjoint(i,l.size());
		r = SpatialDistance(l[i],l[k]);
		t = DX(l,k,i,Dim-1);
		if (mu<(Dim-1)) {
			res += DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DX(l,k,i,mu)*DG2ThermalDrOnr(r,t,beta,a);
		}
		else {
			res += DXDisjoint(l,pk,k,nu,beta)*DXDisjoint(l,pi,i,nu,beta)*DG2ThermalDt(r,t,beta,a);
		}	
	}
	
	v[Dim*j+mu] += f*2.0*res;
}

// PGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWPGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWPGaussianThermalDisjoint_nr(l,j,mu,k,beta,a,f,v);
}

// PGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWPGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWPGaussianThermal2Disjoint_nr(l,j,mu,k,beta,a,f,v);
}

// PRGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWPRGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWPRGaussianThermalDisjoint_nr(l,j,mu,k,nu,beta,a,f,v);
}

// PRGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWPRGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v) {
	uint N = l.size();
	if ((j<N/2 && k<N/2) || (j>=N/2 && k>=N/2))
		NEWPRGaussianThermal2Disjoint_nr(l,j,mu,k,nu,beta,a,f,v);
}

/*----------------------------------------------------------------------------------------------------------------------------
	2 - NEWloopToVector, vectorToLoop
----------------------------------------------------------------------------------------------------------------------------*/

// loopToVector
template<uint Dim>
void NEWloopToVector(const Loop<Dim>& l, vec& v) {
	v.resize(Dim*l.size());
	uint j,k;
	for (j=0; j<l.size(); j++) {
		for (k=0; k<Dim; k++)
			v[j*Dim+k] = (l[j])[k];
	}
}

// vectorToLoop
template<uint Dim>
void NEWvectorToLoop(const vec& v, Loop<Dim>& l) {
	if (v.size()<Dim*l.size()) {
		cerr << "NEWvectorToLoop error: v.size() = " << v.size() << ", loop needs vector of size = " << Dim*l.size() << endl;
	}
	uint j,k;
	for (j=0; j<l.size(); j++) {
		for (k=0; k<Dim; k++)
			(l[j])[k] = v[j*Dim+k];
	}
}

// printAsLoop
void NEWprintAsLoop(const string& f, const uint& Dim, const vec& v) {
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
		cerr << "NEWprintAsLoop error: v.size()%Dim!=0" << endl;
	}
}

// printAsLoop
void NEWprintAsLoop(const string& f, const uint& Dim, const vec& v, const uint len) {
	if (len%Dim==0) {
		uint N = len/Dim, mu;
		ofstream os;
		os.open(f.c_str());
		if (!os.good()) {
			cerr << "NEWprintAsLoop error: " << f << " not opened properly" << endl;
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
		cerr << "NEWprintAsLoop error: len%Dim!=0" << endl;
	}
}

// dimReduce
void NEWdimReduce(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout) {
	if (dimout>dimin) {
		cout << "dimReduce error: dimin=" << dimin << "<dimout=" << dimout << endl;
		return;
	}
	uint NTout = dimout*Nin+zmout, nu;
	vout = Eigen::VectorXd::Zero(NTout);
	for (uint j=0; j<Nin; j++) {
		for (uint mu=0; mu<dimout; mu++) {
			nu = mu+(dimin-dimout);
			vout[j*dimout+mu] = vin[j*dimin+nu];
		}
	}
	for (uint z=0; z<zmout; z++) {
		vout[Nin*dimout+z] = 1.0e-4;
	}
}

// dimIncrease
void NEWdimIncrease(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout) {
	if (dimout<dimin) {
		cout << "dimReduce error: dimin=" << dimin << ">dimout=" << dimout << endl;
		return;
	}
	uint NTout = dimout*Nin+zmout, nu;
	vout = Eigen::VectorXd::Zero(NTout);
	for (uint j=0; j<Nin; j++) {
		for (uint mu=0; mu<dimin; mu++) {
			nu = mu+(dimout-dimin);
			vout[j*dimout+nu] = vin[j*dimin+mu];
		}
	}
	for (uint z=0; z<zmout; z++) {
		vout[Nin*dimout+z] = 1.0e-4;
	}
}


/*----------------------------------------------------------------------------------------------------------------------------
	3 - filename functions
----------------------------------------------------------------------------------------------------------------------------*/

// filenameLoop
template<uint Dim>
Filename NEWfilenameLoopNR(const Parameters& p, const string& base) {
	Filename f = base+"data/nr/loops/dim_"+nts(Dim)+"/K_"+nts(p.K)+"/loop_kappa_"+nts(pow(p.G,3)*p.B)+"_E_"+nts(p.P4)\
		+"_a_"+nts(p.Epsi)+"_mu_"+nts(p.Mu)+".dat";
	if (p.Ng>0)
		(f.Extras).push_back(StringPair("Ng",nts(p.Ng)));
	return f;
}

// filenameThermal
template<uint Dim>
Filename NEWfilenameThermalNR(const Parameters& p, const string& base) {
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
template void NEWL<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWLDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWDistPow<2>(const uint& j, const Loop<2>& l, const number& w, const number& f, number& result);
template void NEWS0<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWS0Disjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWI0<2> (const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWI0Disjoint<2> (const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWIn<2> (const uint& j, const Loop<2>& l, const number& n, const number& f, number& result);
template void NEWInDisjoint<2> (const uint& j, const Loop<2>& l, const number& n, const number& beta, const number& f, number& result);
template void NEWSm<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWAngle<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWFGamma<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWGaussian<2>(const uint& j, const uint& k, const Loop<2>& l, const number& a, const number& f, number& result);
template void NEWGaussianDisjoint<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianLRDisjoint<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermalDisjoint<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermalLRDisjoint<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal2<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal2Disjoint<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal2LRDisjoint<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWMaxXn<2>(const uint& j, const uint& k, const Loop<2>& l, const uint& n, const number& f, number& result);
template void NEWMinXnDisjoint<2>(const uint& j, const uint& k, const Loop<2>& l, const uint& n, const number& f, number& result);
template void NEWInlineCurvatureMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWInlineCurvatureMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWInlineCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWInlineCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWCuspCurvatureMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWCuspCurvatureMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWCuspCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWCuspCurvatureMaxDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWKGMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWKGMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWKGMaxPlane<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWKGMaxPlane<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWKGMaxDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWKGMaxDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWKGMaxPlaneDisjoint<2>(const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWKGMaxPlaneDisjoint<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWmdPX_nr<2>(const Loop<2>& l, const uint& loc, const Point<2>& P, const number& f, vec& v);
template void NEWmdL_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void NEWddL_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void NEWmdLDisjoint_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& beta, const number& p, vec& v);
template void NEWddLDisjoint_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
				const number& beta, const number& p, mat& m);
template void NEWmdDistPow_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& w, const number& p, vec& v);
template void NEWddDistPow_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& w, const number& p, mat& m);
template void NEWmdS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void NEWddS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void NEWmdS0Disjoint_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& beta, const number& p, vec& v);
template void NEWddS0Disjoint_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
			const number& p, mat& m);
template void NEWmdsqrtS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& sqrt2s0, const number& p, vec& v);
template void NEWddsqrtS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
								 const number& sqrt2s0, const number& p, mat& m);
template void NEWmdsqrtS0Disjoint_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& sqrt2s0,\
			 const number& beta, const number& p, vec& v);
template void NEWddsqrtS0Disjoint_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
								 const number& sqrt2s0, const number& beta, const number& p, mat& m);
template void NEWmdI0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& f, vec& v);
template void NEWddI0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<2>& l, const number& f, mat& m);	
template void NEWmdI0Disjoint_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& beta, const number& f, vec& v);
template void NEWddI0Disjoint_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& beta, const number& f, mat& m);
template void NEWmdIn_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& n, const number& f, vec& v);
template void NEWddIn_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<2>& l, const number& n, const number& f, mat& m);	
template void NEWmdInDisjoint_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& n, \
							const number& beta, const number& f, vec& v);
template void NEWddInDisjoint_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& n,\
 							const number& beta, const number& f, mat& m);
template void NEWmdFGamma_nr<2>(const Loop<2>& l, const uint& loc, const number& p, vec& v);
template void NEWddFGamma_nr<2>(const Loop<2>& l, const uint& loc, const number& p, mat& v);
template void NEWErgS0_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void NEWErgS0Disjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void NEWErgL_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void NEWErgLDisjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void NEWErgsqrtS0_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& sqrt2s0, const number& f, number& erg);
template void NEWErgsqrtS0Disjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& sqrt2s0, const number& beta, const number& f, number& erg);
template void NEWErgI0_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void NEWErgI0Disjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void NEWErgIn_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& n, const number& f, number& erg);
template void NEWErgInDisjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, number& erg);
template void NEWErgVor_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg);
template void NEWErgVthr_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgVthrDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgVthrDisjointLR_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussian_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, number& erg);
template void NEWErgGaussianDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianLRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermalDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermalLRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal2_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal2Disjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal2LRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWPS0_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void NEWPS0Disjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void NEWPL_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void NEWPLDisjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void NEWPsqrtS0_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& sqrt2s0, const number& f, vec& v);
template void NEWPsqrtS0Disjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& sqrt2s0, const number& beta, const number& f, vec& v);
template void NEWPI0_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void NEWPI0Disjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void NEWPIn_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& n, const number& f, vec& v);
template void NEWPInDisjoint_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, vec& v);
template void NEWPVor_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v);
template void NEWPRVor_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);
template void NEWPVthr_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRVthr_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPVthrDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRVthrDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussian_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, vec& v);
template void NEWPRGaussian_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);
template void NEWPGaussianDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianLRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianLRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermalDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermalDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermalLRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermalLRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal2_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal2_nr<2>(const Loop<2>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal2Disjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal2Disjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal2LRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal2LRDisjoint_nr<2>(const Loop<2>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWloopToVector<2>(const Loop<2>&,vec&);
template void NEWvectorToLoop<2>(const vec&, Loop<2>&);
template Filename NEWfilenameLoopNR<2>(const Parameters& p, const string&);
template Filename NEWfilenameThermalNR<2>(const Parameters& p, const string&);
template void NEWPseudoAngle<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void NEWVor<2>(const uint& j, const uint& k, const Loop<2>& l, const number& a, const number& f, number& result);
template void NEWVlr<2>(const uint& j, const uint& k, const Loop<2>& l, const number& a, const number& f, number& result);
template void NEWVer<2>(const uint& j, const uint& k, const Loop<2>& l, const number& a, const number& f, number& result);
template void NEWVdr<2>(const uint& j, const uint& k, const Loop<2>& l, const number& a, const number& f, number& result);
template void NEWVthr<2>(const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWVthrDisjoint<2> (const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWVthrDisjointLR<2> (const uint& j, const uint& k, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWVnonrelDisjoint<2> (const uint& j, const Loop<2>& l, const number& beta, const number& f, number& result);
template void NEWVnonrelrDisjoint<2> (const uint& j, const Loop<2>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWmdVor_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, const number& a, const number& f, vec& v);
template void NEWmdVer_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, const number& a, const number& f, vec& v);
template void NEWmdVdr_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, const number& a, const number& f, vec& v);
template void NEWmdVthr_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, \
			const number& beta, const number& a, const number& f, vec& v);
template void NEWmdVnonrelDisjoint_nr<2> (const uint& j, const uint& mu, const Loop<2>& l, const number& beta, const number& f, vec& v);
template void NEWmdVnonrelrDisjoint_nr<2> (const uint& j, const uint& mu, const Loop<2>& l, const number& beta, const number& a, const number& f, vec& v);
template void NEWmdVthrDisjoint_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, \
			const number& beta, const number& a, const number& f, vec& v);
template void NEWmdVthrDisjointLR_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, \
			const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussian_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, \
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal2_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, \
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermalDisjoint_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal2Disjoint_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermalLRDisjoint_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal2LRDisjoint_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWddVor_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& a, const number& f, mat& m);
template void NEWddVer_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& a, const number& f, mat& m);
template void NEWddVdr_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& a, const number& f, mat& m);
template void NEWddVthr_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddVthrDisjoint_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddVthrDisjointLR_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddVnonrelDisjoint_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& f, mat& m);
template void NEWddVnonrelrDisjoint_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal2_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermalDisjoint_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal2Disjoint_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermalLRDisjoint_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal2LRDisjoint_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussian_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& a, const number& f, mat& m);
template void NEWmdGaussianDisjoint_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianLRDisjoint_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWddGaussianDisjoint_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianLRDisjoint_nr<2> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWmdVlr_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, const number& a, const number& f, vec& v);
template void NEWddVlr_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
						 const number& a, const number& f, mat& m);					 
						 



// dim 4
template void NEWL<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWLDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWDistPow<4>(const uint& j, const Loop<4>& l, const number& w, const number& f, number& result);
template void NEWS0<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWS0Disjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWI0<4> (const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWI0Disjoint<4> (const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWIn<4> (const uint& j, const Loop<4>& l, const number& n, const number& f, number& result);
template void NEWInDisjoint<4> (const uint& j, const Loop<4>& l, const number& n, const number& beta, const number& f, number& result);
template void NEWSm<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWAngle<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWFGamma<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWGaussian<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result);
template void NEWGaussianDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianLRDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermalDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermalLRDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal2<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal2Disjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWGaussianThermal2LRDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWMaxXn<4>(const uint& j, const uint& k, const Loop<4>& l, const uint& n, const number& f, number& result);
template void NEWMinXnDisjoint<4>(const uint& j, const uint& k, const Loop<4>& l, const uint& n, const number& f, number& result);
template void NEWInlineCurvatureMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWInlineCurvatureMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWInlineCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWInlineCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWCuspCurvatureMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWCuspCurvatureMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWCuspCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWCuspCurvatureMaxDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWKGMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWKGMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWKGMaxPlane<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWKGMaxPlane<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void NEWKGMaxDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWKGMaxDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWKGMaxPlaneDisjoint<4>(const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWKGMaxPlaneDisjoint<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);
template void NEWmdPX_nr<4>(const Loop<4>& l, const uint& loc, const Point<4>& P, const number& f, vec& v);
template void NEWmdL_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void NEWddL_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void NEWmdLDisjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& p, vec& v);
template void NEWddLDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
				const number& beta, const number& p, mat& m);
template void NEWmdDistPow_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& w, const number& p, vec& v);
template void NEWddDistPow_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& w, const number& p, mat& m);
template void NEWmdS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void NEWddS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void NEWmdS0Disjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& p, vec& v);
template void NEWddS0Disjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
			const number& p, mat& m);
template void NEWmdsqrtS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& sqrt4s0, const number& p, vec& v);
template void NEWddsqrtS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
								 const number& sqrt4s0, const number& p, mat& m);
template void NEWmdsqrtS0Disjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& sqrt4s0,\
			 const number& beta, const number& p, vec& v);
template void NEWddsqrtS0Disjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
								 const number& sqrt4s0, const number& beta, const number& p, mat& m);
template void NEWmdI0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& f, vec& v);
template void NEWddI0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<4>& l, const number& f, mat& m);	
template void NEWmdI0Disjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& f, vec& v);
template void NEWddI0Disjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& beta, const number& f, mat& m);
template void NEWmdIn_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& n, const number& f, vec& v);
template void NEWddIn_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<4>& l, const number& n, const number& f, mat& m);	
template void NEWmdInDisjoint_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& n, \
							const number& beta, const number& f, vec& v);
template void NEWddInDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& n,\
 							const number& beta, const number& f, mat& m);
template void NEWmdFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, vec& v);
template void NEWddFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, mat& v);
template void NEWErgS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void NEWErgS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void NEWErgL_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void NEWErgLDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void NEWErgsqrtS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, number& erg);
template void NEWErgsqrtS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, number& erg);
template void NEWErgI0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, number& erg);
template void NEWErgI0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);
template void NEWErgIn_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& f, number& erg);
template void NEWErgInDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, number& erg);
template void NEWErgVor_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg);
template void NEWErgVthr_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgVthrDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgVthrDisjointLR_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussian_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, number& erg);
template void NEWErgGaussianDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermalDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermalLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal2_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal2Disjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWErgGaussianThermal2LRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);
template void NEWPS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void NEWPS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void NEWPL_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void NEWPLDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void NEWPsqrtS0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, vec& v);
template void NEWPsqrtS0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, vec& v);
template void NEWPI0_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& f, vec& v);
template void NEWPI0Disjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);
template void NEWPIn_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& f, vec& v);
template void NEWPInDisjoint_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, vec& v);
template void NEWPVor_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v);
template void NEWPRVor_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);
template void NEWPVthr_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRVthr_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPVthrDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRVthrDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussian_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, vec& v);
template void NEWPRGaussian_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);
template void NEWPGaussianDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermalDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermalDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermalLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermalLRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal2_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal2_nr<4>(const Loop<4>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal2Disjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal2Disjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWPGaussianThermal2LRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);
template void NEWPRGaussianThermal2LRDisjoint_nr<4>(const Loop<4>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);
template void NEWloopToVector<4>(const Loop<4>&,vec&);
template void NEWvectorToLoop<4>(const vec&, Loop<4>&);
template Filename NEWfilenameLoopNR<4>(const Parameters& p, const string&);
template Filename NEWfilenameThermalNR<4>(const Parameters& p, const string&);
template void NEWPseudoAngle<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void NEWVor<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result);
template void NEWVlr<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result);
template void NEWVer<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result);
template void NEWVdr<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result);
template void NEWVthr<4>(const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWVthrDisjoint<4> (const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWVthrDisjointLR<4> (const uint& j, const uint& k, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWVnonrelDisjoint<4> (const uint& j, const Loop<4>& l, const number& beta, const number& f, number& result);
template void NEWVnonrelrDisjoint<4> (const uint& j, const Loop<4>& l, const number& beta, const number& a, const number& f, number& result);
template void NEWmdVor_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v);
template void NEWmdVer_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v);
template void NEWmdVdr_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v);
template void NEWmdVthr_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
			const number& beta, const number& a, const number& f, vec& v);
template void NEWmdVnonrelDisjoint_nr<4> (const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& f, vec& v);
template void NEWmdVnonrelrDisjoint_nr<4> (const uint& j, const uint& mu, const Loop<4>& l, const number& beta, const number& a, const number& f, vec& v);
template void NEWmdVthrDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
			const number& beta, const number& a, const number& f, vec& v);
template void NEWmdVthrDisjointLR_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
			const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussian_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal2_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, \
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermalDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal2Disjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermalLRDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianThermal2LRDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWddVor_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m);
template void NEWddVer_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m);
template void NEWddVdr_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m);
template void NEWddVthr_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddVthrDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddVthrDisjointLR_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddVnonrelDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& f, mat& m);
template void NEWddVnonrelrDisjoint_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal2_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermalDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal2Disjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermalLRDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianThermal2LRDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussian_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m);
template void NEWmdGaussianDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWmdGaussianLRDisjoint_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l,\
				const number& beta, const number& a, const number& f, vec& v);
template void NEWddGaussianDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWddGaussianLRDisjoint_nr<4> (const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
template void NEWmdVlr_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v);
template void NEWddVlr_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
						 const number& a, const number& f, mat& m);
