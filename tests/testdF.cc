/*
	testdF
		writing a function to take as an argument the scalar functions dFdt etc and then to spit our mdF as a vector
*/

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "analysis.h"
#include "check.h"
#include "folder.h"
#include "genloop.h"
#include "nrloop.h"
#include "print.h"
#include "simple.h"

using namespace std;

typedef number(*ThermalFunction)(const number&, const number&, const number&, const number&);

typedef number(*VacuumFunction)(const number&, const number&);

// F0
number F0(const number& x, const number& a) {
	return  -1.0/(4.0*pow(PI,2)*(pow(a,2) + pow(x,2)));
}

// DF0DrOnr
number DF0DrOnr(const number& x, const number& a) {
	return 1.0/(2.0*pow(PI,2)*pow(pow(a,2) + pow(x,2),2));
}
 
 // DDF0DrDr
number DDF0DrDr(const number& x, const number& a) {
	return (pow(a,2) - 3.0*pow(x,2))/(2.0*pow(PI,2)*pow(pow(a,2) + pow(x,2),3));
}

// V0Generic
template <uint Dim>
static void V0Generic (const uint& j, const uint& k, const Loop<Dim>& l, VacuumFunction FVacuum, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
		number x_jk = Distance(l[j],l[k]);
		
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*FVacuum(x_jk,a);	
	}
	
}
 
// mdF0r_nr
template<uint Dim>
void mdF0_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
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
		res += (-1.0/pow(2.0*PI,2))*2.0*(l[j])[mu]/a/a; 
	if (i==mj)
		res += (-1.0/pow(2.0*PI,2))*(-(l[mj])[mu]/a/a);
	if (i==pj)
		res += (-1.0/pow(2.0*PI,2))*(-(l[pj])[mu]/a/a);
		
	v[j*Dim+mu] += -f*(-pow(2.0*PI,2))*res;
}

// ddV0r_nr
template<uint Dim>
void ddV0r_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
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
		res += (-1.0/pow(2.0*PI,2))*(2.0/a/a);
	if (k==mj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	if (k==pj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	
	m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
	
}

// FThermal
number FThermal(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return  sinh((2.0*PI*ra)/beta)/(4.0*PI*ra*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))); //
}

// DFThermalDrOnr
number DFThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return cosh((2.0*PI*ra)/beta)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 - sinh((2.0*PI*ra)/beta)/\
 (4.0*PI*pow(ra,3)*beta*(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta))) \
 + pow(sinh((2.0*PI*ra)/beta),2)/\
 (2.0*(pow(a,2) + pow(r,2))*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2));
 }

// DFThermalDt
number DFThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (sin((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (2.0*ra*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2));
}
 

// DDFThermalDrDr
number DDFThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
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
number DDFThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return (PI*cosh((2.0*PI*ra)/beta)*sin((2.0*PI*t)/beta))/\
 ((pow(a,2) + pow(r,2))*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 - (sin((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (2.0*pow(ra,2*1.5)*pow(beta,2)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2))\
 + (2.0*PI*sin((2.0*PI*t)/beta)*pow(sinh((2.0*PI*ra)/beta),2))/\
 ((pow(a,2) + pow(r,2))*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),3));
}

// DDFThermalDtDt
number DDFThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	number ra = sqrt(pow(a,2) + pow(r,2));
	return  (PI*cos((2.0*PI*t)/beta)*sinh((2.0*PI*ra)/beta))/\
 (ra*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),2)) \
 + (2.0*PI*pow(sin((2.0*PI*t)/beta),2)*sinh((2.0*PI*ra)/beta))/\
 (ra*pow(beta,3)*pow(cos((2.0*PI*t)/beta) - cosh((2.0*PI*ra)/beta),3));
}

// VthGeneric
template <uint Dim>
static void VthGeneric (const uint& j, const uint& k, const Loop<Dim>& l, ThermalFunction Fth, const number& beta, const number& a, const number& f, number& result) {

	if (k<=j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
		number r = SpatialDistance(l[j],l[k]);
		number t = DX(l,k,j,Dim-1);
		
		result += f*(-pow(2.0*PI,2))*(1.0+(number)(k<j))*Dot(l[pj],l[j],l[pk],l[k])*Fth(r,t,beta,a);	
	}
	
}

 // mdVthr_nr
template<uint Dim>
void mdFth_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					ThermalFunction Fth, ThermalFunction dFthdronr, ThermalFunction dFthdt, \
					 	const number& beta, const number& a, const number& f, vec& v) {
	number res = 0.0;
	
	uint pj = posNeigh(j,l.size());
	uint mj = negNeigh(j,l.size());
	uint pi = posNeigh(i,l.size());
		
	if (i!=j) {
		number r_ij = SpatialDistance(l[i],l[j]);
		number t_ij = DX(l,j,i,Dim-1); // checked order
		number FThermal_ij = Fth(r_ij,t_ij,beta,a);
		number DFThermalDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
		number DFThermalDt_ij = dFthdt(r_ij,t_ij,beta,a);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		
		res += 2.0*FThermal_ij*(-DX(l,i,mu));
		if (mu<(Dim-1))
			res += 2.0*DFThermalDrOnr_ij*DX(l,j,i,mu)*T_ij;
		else
			res += 2.0*DFThermalDt_ij*T_ij;
	}
	
	if (i!=mj) {
		number r_imj = SpatialDistance(l[i],l[mj]);
		number t_imj = DX(l,mj,i,Dim-1); // checked order
		number FThermal_imj = Fth(r_imj,t_imj,beta,a);
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
		
	v[j*Dim+mu] += -f*(-pow(2.0*PI,2))*res;
}

// ddVthr_nr
template <uint Dim>
void ddVthr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
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
			number DFThermalDrOnr_jk = dFthdronr(r_jk,t_jk,beta,a);
			number DFThermalDt_jk = dFthdt(r_jk,t_jk,beta,a);
			number DDFThermalDrDr_jk = ddFthdrdr(r_jk,t_jk,beta,a);
			number DDFThermalDrDtOnr_jk = ddFthdrdtonr(r_jk,t_jk,beta,a);
			number DDFThermalDtDt_jk = ddFthdtdt(r_jk,t_jk,beta,a);
			if (mu==(Dim-1)) {
				res += - 2.0*DFThermalDt_jk*(-DX(l,j,nu)); //
				if (nu<(Dim-1))
					res += (2.0*DDFThermalDrDtOnr_jk*DX(l,j,k,nu)*T_jk); //
				else
					res += -2.0*DDFThermalDtDt_jk*T_jk;//
			}
			else {
				res += + 2.0*DFThermalDrOnr_jk*(-DX(l,j,nu))*DX(l,j,k,mu); //
				if (nu==(Dim-1)) 
					res += + 2.0*DDFThermalDrDtOnr_jk*DX(l,j,k,mu)*T_jk; //
				else if (r_jk>MIN_NUMBER) //as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
					res += (2.0*DFThermalDrOnr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2) \
							- (2.0*DDFThermalDrDr_jk*DX(l,j,k,mu)*DX(l,j,k,nu)*T_jk)/pow(r_jk,2); //
				if (nu==mu)
					res += - (2.0*DFThermalDrOnr_jk*T_jk);//
			}
			
			if (nu==(Dim-1)) {
				res += 2.0*DFThermalDt_jk*(-DX(l,k,mu)); //
			}
			else {
				res += - 2.0*DFThermalDrOnr_jk*DX(l,j,k,nu)*(-DX(l,k,mu)); //
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
			t_ij = DX(l,j,i,Dim-1);
			t_imj = DX(l,mj,i,Dim-1);

			if (k==j && i!=j) {
				DFThermalDrOnr_ij = dFthdronr(r_ij,t_ij,beta,a);
				DFThermalDt_ij = dFthdt(r_ij,t_ij,beta,a);
				DDFThermalDrDr_ij = ddFthdrdr(r_ij,t_ij,beta,a);
				DDFThermalDrDtOnr_ij = ddFthdrdtonr(r_ij,t_ij,beta,a);
				DDFThermalDtDt_ij = ddFthdtdt(r_ij,t_ij,beta,a);
				if (mu==(Dim-1)) {
					res += 2.0*DFThermalDt_ij*(-DX(l,i,nu)); //
					if (nu==(Dim-1)) 
						res += 2.0*DDFThermalDtDt_ij*T_ij; //
					else
						res += 2.0*DDFThermalDrDtOnr_ij*DX(l,j,i,nu)*T_ij; //
				}
				else {
					res += 2.0*DFThermalDrOnr_ij*(-DX(l,i,nu))*DX(l,j,i,mu); //
					if (nu==(Dim-1))
						res += 2.0*DDFThermalDrDtOnr_ij*DX(l,j,i,mu)*T_ij; //
					else if (r_ij>MIN_NUMBER) { // as -DFThermalDrOnr+DDFThermalDrDr->0 as r->0
						res += - (2.0*DFThermalDrOnr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2) \
							+ (2.0*DDFThermalDrDr_ij*DX(l,j,i,mu)*DX(l,j,i,nu)*T_ij)/pow(r_ij,2);  //
					}
					if (nu==mu)
						res += 2.0*DFThermalDrOnr_ij*T_ij; //
				}
				if (nu==(Dim-1)) {
					res += 2.0*DFThermalDt_ij*(-DX(l,i,mu)); //
				}
				else {
					res += - 2.0*DFThermalDrOnr_ij*(-DX(l,i,mu))*DX(l,i,j,nu); //
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
		res += (-1.0/pow(2.0*PI,2))*2.0/a/a;
	if (k==mj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	if (k==pj && mu==nu)
		res += (-1.0/pow(2.0*PI,2))*(-1.0/a/a);
	
	m(Dim*j+mu,Dim*k+nu) += f*(-pow(2.0*PI,2))*res;
}

void dimReduce(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout) {
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

void dimIncrease(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout) {
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

int main() {

cout << "testdF:" << endl;

// inputs
string inputsFile = "inputs4";
cout << "loading parameters from " << inputsFile << endl;
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min;
if (p.empty()) {
	cerr << "Parameters empty: nothing in inputs file: " << inputsFile << endl;
	return 1;
}

#define dim 4
#define dimred 2
uint N = pow(2,p.K);
uint zm = dim+1;
uint zmred = dimred+1;
uint NT = N*dim+zm;
uint NTred = N*dimred+zmred;

// defining vector and matrix quantities
vec x(N*dim), xred(N*dimred);
vec mds(NT), mdsred(NTred);
mat dds(NT,NT);

// defining xLoop
Loop<dim> xLoop(p.K,0);
Loop<dimred> xLoopRed(p.K,0);

string loadFile = "data/nr/loops/dim_4/K_6/loop_kappa_1_E_0_a_0.1_mu_0.1.dat";
cout << "loading x from " << loadFile << endl;

// loading x
loadVectorBinary(loadFile,x);
if (x.size()<NT) {
	x.conservativeResize(NT);
	for (uint mu=0; mu<zm; mu++)
		x[N*dim+mu] = 1.0e-4;
}
else if (x.size()>NT)
	x.conservativeResize(NT);
	
// loading x to xLoop
vectorToLoop(x,xLoop);

// dimReduce and dimIncrease test
dimReduce(x,dim,N,xred,dimred,zmred);
vec xTest(N*dim);
dimIncrease(xred,dimred,N,xTest,dim,zm);
cout << "dimReduce then dimIncrease:" << endl;
cout << "norm before:" << x.norm() << endl;
cout << "norm reduced:" << xred.norm() << endl;
cout << "norm after:" << xTest.norm() << endl;

vectorToLoop(xred,xLoopRed);

number beta = 1.0/p.T;
number g = p.G*p.G/8.0/PI;
number v = 0.0;

clock_t time;

cout << "-------------------------------------------------------------------------------------------------------------------------" << endl;
cout << "dim = "<< dim << endl;

cout << "v before = " << v << endl;
cout << "norm(mds) before: " << mds.norm() << endl;
cout << "norm(dds) before: " << dds.norm() << endl;

time = clock();

// bulk
for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dim; mu++) {
		for (uint k=0; k<N; k++) {
			if (mu==0) {
				V0Generic(j,k,xLoop,&F0,p.Epsi,g,v);
				VthGeneric(j,k,xLoop,&FThermal,beta,p.Epsi,g,v);
			}
			mdF0_nr(j, mu, k, xLoop, &F0, &DF0DrOnr, p.Epsi, g, mds);
			mdFth_nr(j, mu, k, xLoop, &FThermal, &DFThermalDrOnr, &DFThermalDt, beta, p.Epsi, g, mds);
			for (uint nu=0; nu<dim; nu++) {
				ddV0r_nr(j, mu, k, nu, xLoop, &F0, &DF0DrOnr, &DDF0DrDr, p.Epsi, g, dds);
				ddVthr_nr(j, mu, k, nu, xLoop, &FThermal, &DFThermalDrOnr, &DFThermalDt, &DDFThermalDrDr, &DDFThermalDtDrOnr, &DDFThermalDtDt, beta, p.Epsi, g, dds);
			}
		}
	}
}

cout << "v middle = " << v << endl;
cout << "norm(mds) middle: " << mds.norm() << endl;
cout << "norm(dds) middle: " << dds.norm() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;
time = clock();

// bulk
for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dim; mu++) {
		for (uint k=0; k<N; k++) {
			if (mu==0) {
				Vor(j,k,xLoop,p.Epsi,-g,v);
				Vthr(j,k,xLoop,beta,p.Epsi,-g,v);
			}
			mdVor_nr(j, mu, k, xLoop, p.Epsi, -g, mds);
			mdVthr_nr(j, mu, k, xLoop, beta, p.Epsi, -g, mds);
			for (uint nu=0; nu<dim; nu++) {
				ddVor_nr(j, mu, k, nu, xLoop, p.Epsi, -g, dds);
				ddVthr_nr(j, mu, k, nu, xLoop, beta, p.Epsi, -g, dds);
			}
		}
	}
}

cout << "v after = " << v << endl;
cout << "norm(mds) after: " << mds.norm() << endl;
cout << "norm(dds) after: " << dds.norm() << endl;
cout << "mds.size() = " << mds.size() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;

cout << "before and after should be zero, middle nonzero." << endl << endl;

cout << "-------------------------------------------------------------------------------------------------------------------------" << endl;
cout << "dimred = "<< dimred << endl;

cout << "norm(mdsred) before: " << mdsred.norm() << endl;

time = clock();

// bulk
for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dimred; mu++) {
		for (uint k=0; k<N; k++) {
			mdF0_nr(j, mu, k, xLoopRed, &F0, &DF0DrOnr, p.Epsi, g, mdsred);
			mdFth_nr(j, mu, k, xLoopRed, &FThermal, &DFThermalDrOnr, &DFThermalDt, beta, p.Epsi, g, mdsred);
		}
	}
}

cout << "norm(mdsred) middle: " << mdsred.norm() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;
time = clock();

// bulk
/*for (uint j=0; j<N; j++) {
	for (uint mu=0; mu<dimred; mu++) {
		for (uint k=0; k<N; k++) {
			mdVthr_nr(j, mu, k, xLoopRed, beta, p.Epsi, -g, mdsred);
		}
	}
}*/
cout << "dropped this bit as don't have the code yet" << endl;

cout << "norm(mdsred) after: " << mdsred.norm() << endl;
cout << "mdsred.size() = " << mdsred.size() << endl;

//stopping clock
time = clock() - time;
cout << "time taken: "  << time/1000000.0 << endl;
cout << "-------------------------------------------------------------------------------------------------------------------------" << endl;

return 0;
}
