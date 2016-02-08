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
		1 - nr loop functions
		2 - loopToVector, vectorToLoop
		3 - explicit instatiation
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - nr loop functions
----------------------------------------------------------------------------------------------------------------------------*/

// L
template <uint Dim>
void L (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pos = (j==(l.size()-1)? 0: j+1);
	result += f*Distance(l[pos],l[j]);
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
	result += f*DistanceSquared(l[pj],l[j])/4.0;
}

// I0
template <uint Dim>
void I0 (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	cerr << "mdI_nr error: not defined for dim = " << Dim << endl;
}

// Gamma
template <uint Dim>
void Gamma (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	number cos_gamma;
	
	cos_gamma = Dot(l[pj],l[j],l[j],l[nj]);
	cos_gamma /= Distance(l[pj],l[j])*Distance(l[j],l[nj]);
	
	result += f*(acos(cos_gamma));
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

// V1r
template <uint Dim>
void V1r (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])*pow(DistanceSquared(l[j],l[k])+a*a,(2.0-Dim)/2.0);	
	}
}

// V2r
template <uint Dim>
void V2r (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])*pow(0.25*DistanceSquared(l[pj]+l[j],l[pk]+l[k])+a*a,(2.0-Dim)/2.0);	
	}
}

// SimpleCurvatureMax
template <uint Dim>
void SimpleCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1):j-1);
	Point<Dim> dd = l[pj] + l[nj] - 2.0*l[j];
	Point<Dim> d = l[pj] - l[j];
	number temp = f*abs(Dot(dd,d))/Dot(d,d);
	result = (temp>result? temp: result);
}

// SimpleCurvatureMax
template <uint Dim>
void SimpleCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint nj = (j==0? (l.size()-1):j-1);
		Point<Dim> dd = l[pj] + l[nj] - 2.0*l[j];
		Point<Dim> d = l[pj] - l[j];
		number temp = f*abs(Dot(dd,d))/Dot(d,d);
		result = (temp>result? temp: result);
	}
}

// SimpleCurvatureAvg
template <uint Dim>
void SimpleCurvatureAvg (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1):j-1);
	Point<Dim> dd = l[pj] + l[nj] - 2.0*l[j];
	Point<Dim> d = l[pj] - l[j];
	number temp = f*abs(Dot(dd,d))/Dot(d,d);
	result += temp/(number)l.size();
}

// SimpleCurvatureAvg
template <uint Dim>
void SimpleCurvatureAvg (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint nj = (j==0? (l.size()-1):j-1);
		Point<Dim> dd = l[pj] + l[nj] - 2.0*l[j];
		Point<Dim> d = l[pj] - l[j];
		number temp = f*abs(Dot(dd,d))/Dot(d,d);
		result += temp/((number)l.size()-2.0);
	}
}

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pos1, pos2;
	Point<Dim> pp;
	number t_sqrd, t_dot_sqrd, t_t_dot, temp;
	pos1 = (j==(l.size()-1)? 0 : j+1);
	pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
	pp = l[pos1];
	pp *= 2.0;
	pp -= l[pos2];
	t_sqrd = DistanceSquared(l[pos1],l[j]);
	t_dot_sqrd = DistanceSquared(l[j],pp);
	t_t_dot = Dot(l[j],pp,l[pos1],l[j]);
	temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
	result =  (temp>result? temp: result);
}

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pos1, pos2;
		Point<Dim> pp;
		number t_sqrd, t_dot_sqrd, t_t_dot, temp;
		pos1 = (j==(l.size()-1)? 0 : j+1);
		pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
		pp = l[pos1];
		pp *= 2.0;
		pp -= l[pos2];
		t_sqrd = DistanceSquared(l[pos1],l[j]);
		t_dot_sqrd = DistanceSquared(l[j],pp);
		t_t_dot = Dot(l[j],pp,l[pos1],l[j]);
		temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size();
		result = (temp>result? temp: result);
	}
}

// KGAvg
template <uint Dim>
void KGAvg (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pos1, pos2;
	Point<Dim> pp;
	number t_sqrd, t_dot_sqrd, t_t_dot, temp;
	pos1 = (j==(l.size()-1)? 0 : j+1);
	pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
	pp = l[pos1];
	pp *= 2.0;
	pp -= l[pos2];
	t_sqrd = DistanceSquared(l[pos1],l[j]);
	t_dot_sqrd = DistanceSquared(l[j],pp);
	t_t_dot = Dot(l[j],pp,l[pos1],l[j]);
	temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd);
	result +=  temp;
}

// KGAvg
template <uint Dim>
void KGAvg (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pos1, pos2;
		Point<Dim> pp;
		number t_sqrd, t_dot_sqrd, t_t_dot, temp;
		pos1 = (j==(l.size()-1)? 0 : j+1);
		pos2 = (pos1==(l.size()-1)? 0 : pos1+1);
		pp = l[pos1];
		pp *= 2.0;
		pp -= l[pos2];
		t_sqrd = DistanceSquared(l[pos1],l[j]);
		t_dot_sqrd = DistanceSquared(l[j],pp);
		t_t_dot = Dot(l[j],pp,l[pos1],l[j]);
		temp = f*sqrt((t_dot_sqrd-t_t_dot*t_t_dot/t_sqrd)/t_sqrd)*l.size()/((number)l.size()-2.0);
		result += temp;
	}
}

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1):j-1);
	number dx = (l[pj])[2]-(l[j])[2];
	number ddx = (l[pj])[2]-2.0*(l[j])[2]+(l[nj])[2];
	number dy = (l[pj])[3]-(l[j])[3];
	number ddy = (l[pj])[3]-2.0*(l[j])[3]+(l[nj])[3];
	number temp = f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5);
	result = (temp>result? temp: result);
}

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint nj = (j==0? (l.size()-1):j-1);
		number dx = (l[pj])[2]-(l[j])[2];
		number ddx = (l[pj])[2]-2.0*(l[j])[2]+(l[nj])[2];
		number dy = (l[pj])[3]-(l[j])[3];
		number ddy = (l[pj])[3]-2.0*(l[j])[3]+(l[nj])[3];
		number temp = f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5);
		result = (temp>result? temp: result);
	}
}

// KGAvgPlane
template <uint Dim>
void KGAvgPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1):j-1);
	number dx = (l[pj])[2]-(l[j])[2];
	number ddx = (l[pj])[2]-2.0*(l[j])[2]+(l[nj])[2];
	number dy = (l[pj])[3]-(l[j])[3];
	number ddy = (l[pj])[3]-2.0*(l[j])[3]+(l[nj])[3];
	result += f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5)/(number)l.size();
}

// KGAvgPlane
template <uint Dim>
void KGAvgPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result) {
	if (j!=ex1 && j!=ex2) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint nj = (j==0? (l.size()-1):j-1);
		number dx = (l[pj])[2]-(l[j])[2];
		number ddx = (l[pj])[2]-2.0*(l[j])[2]+(l[nj])[2];
		number dy = (l[pj])[3]-(l[j])[3];
		number ddy = (l[pj])[3]-2.0*(l[j])[3]+(l[nj])[3];
		result += f*abs(dx*ddy-dy*ddx)/pow(dx*dx+dy*dy,1.5)/((number)l.size()-2.0);
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
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	
	number norm = Distance(l[j],l[nj]);
	number temp = ((l[j])[mu]-(l[nj])[mu])/norm;
	
	norm = Distance(l[pj],l[j]);
	temp += ((l[j])[mu]-(l[pj])[mu])/norm;
	
	v[j*Dim+mu] += -f*temp;
}

// ddL_nr
template<uint Dim>
void ddL_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	if (k==j) {
		number normn = Distance(l[j],l[nj]), normp = Distance(l[pj],l[j]), temp = 0.0;
		if (mu==nu)
			temp += 1.0/normn + 1.0/normp;
		temp -= ((l[j])[mu]-(l[nj])[mu])*((l[j])[nu]-(l[nj])[nu])/pow(normn,3);
		temp -= ((l[j])[mu]-(l[pj])[mu])*((l[j])[nu]-(l[pj])[nu])/pow(normp,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==nj) {
		number norm = Distance(l[j],l[nj]), temp = 0.0;
		if (mu==nu)
			temp -= 1.0/norm;
		temp += ((l[j])[mu]-(l[nj])[mu])*((l[j])[nu]-(l[nj])[nu])/pow(norm,3);
		m(j*Dim+mu,k*Dim+nu) += f*temp;
	}
	else if (k==pj) {
		number norm = Distance(l[pj],l[j]), temp = 0.0;
		if (mu==nu)
			temp -= 1.0/norm;
		temp += ((l[j])[mu]-(l[pj])[mu])*((l[j])[nu]-(l[pj])[nu])/pow(norm,3);
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

// mdS0_nr
template<uint Dim>
void mdS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	v[j*Dim+mu] += -f*(2.0*(l[j])[mu] - (l[nj])[mu] - (l[pj])[mu])*(number)l.size()/2.0;
}

// ddS0_nr
template<uint Dim>
void ddS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m) {
	if (mu==nu) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint nj = (j==0? (l.size()-1): j-1);
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
	
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	uint pk = (k==(l.size()-1)? 0: k+1);
	uint nk = (k==0? (l.size()-1): k-1);
	
	m(j*Dim+mu,k*Dim+nu) -= (f*pow((number)l.size(),2)/pow(sqrt4s0,3)) * (2.0*(l[j])[mu]-(l[pj])[mu]-(l[nj])[mu]) \
												* (2.0*(l[k])[nu]-(l[pk])[nu]-(l[nk])[nu]);
	
}

// mdV1r_nr
template<uint Dim>
void mdV1r_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdV1r_nr( j, mu, i, l, a, f, v);
	
}

// mdV1r_nr
template<uint Dim>
void mdV1r_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	cerr << "mdV1r_nr Error: no script written for dim = " << Dim << endl;
}

// ddV1r_nr
template<uint Dim>
void ddV1r_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "mddV1r_nr Error: no script written for dim = " << Dim << endl;
}

// mdV2r_nr
template<uint Dim>
void mdV2r_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v) {

	for (uint i=0; i<l.size(); i++)
		mdV2r_nr( j, mu, i, l, a, f, v);
	
}

// mdV2r_nr
template<uint Dim>
void mdV2r_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	cerr << "mdV2r_nr Error: no script written for dim = " << Dim << endl;
}

// ddV2r_nr
template<uint Dim>
void ddV2r_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m) {
	cerr << "mddV2r_nr Error: no script written for dim = " << Dim << endl;
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
	if (v.size()!=Dim*l.size() && v.size()!=Dim*(l.size()+1)) {
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
	3 - explicit instantiation
----------------------------------------------------------------------------------------------------------------------------*/

// dim 2
template void L<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void S0<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void Sm<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void Gamma<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void FGamma<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void SimpleCurvatureMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void SimpleCurvatureMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void SimpleCurvatureAvg<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void SimpleCurvatureAvg<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMax<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void KGMax<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGAvg<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void KGAvg<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMaxPlane<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void KGMaxPlane<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGAvgPlane<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void KGAvgPlane<2>(const uint& j, const Loop<2>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void mdPX_nr<2>(const Loop<2>& l, const uint& loc, const Point<2>& P, const number& f, vec& v);
template void mdL_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void ddL_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void mdS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void ddS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void mdsqrtS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& sqrt4s0, const number& p, vec& v);
template void ddsqrtS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
								 const number& sqrt4s0, const number& p, mat& m);
template void mdV1r_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& a, const number& p, vec& v);
template void mdFGamma_nr<2>(const Loop<2>& l, const uint& loc, const number& f, vec& v);
template void ddFGamma_nr<2>(const Loop<2>& l, const uint& loc, const number& f, mat& v);
template void loopToVector<2>(const Loop<2>&,vec&);
template void vectorToLoop<2>(const vec&, Loop<2>&);

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
template void S0<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void Sm<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void Gamma<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void FGamma<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void SimpleCurvatureMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void SimpleCurvatureMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void SimpleCurvatureAvg<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void SimpleCurvatureAvg<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMax<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void KGMax<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGAvg<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void KGAvg<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGMaxPlane<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void KGMaxPlane<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void KGAvgPlane<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void KGAvgPlane<4>(const uint& j, const Loop<4>& l, const uint& ex1, const uint& ex2, const number& f, number& result);
template void mdPX_nr<4>(const Loop<4>& l, const uint& loc, const Point<4>& P, const number& f, vec& v);
template void mdL_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void ddL_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void mdS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void ddS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void mdsqrtS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& sqrt4s0, const number& p, vec& v);
template void ddsqrtS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
								 const number& sqrt4s0, const number& p, mat& m);
template void mdV1r_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& a, const number& p, vec& v);
template void mdFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, vec& v);
template void ddFGamma_nr<4>(const Loop<4>& l, const uint& loc, const number& p, mat& v);
template void loopToVector<4>(const Loop<4>&,vec&);
template void vectorToLoop<4>(const vec&, Loop<4>&);

// V1r
template <> void V1r<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])/(DistanceSquared(l[j],l[k])+a*a);	
	}
}

// V2r
template <> void V2r<4>(const uint& j, const uint& k, const Loop<4>& l, const number& a, const number& f, number& result) {

	if (k<j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint pk = (k==(l.size()-1)? 0: k+1);
	
		result += f*2.0*Dot(l[pj],l[j],l[pk],l[k])/(0.25*DistanceSquared(l[pj]+l[j],l[pk]+l[k])+a*a);	
	}
}

// mdV1r_nr
template <> void mdV1r_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v) {
	number res = 0.0;
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint mj = (j==0? (l.size()-1): j-1);
	uint pi = (i==(l.size()-1)? 0: i+1);
		
	if (i!=j) {
		number B_ij = a*a + DistanceSquared(l[i],l[j]);
		number T_ij = Dot(l[pi],l[i],l[pj],l[j]);
		res += 2.0*DX(l,i,pi,mu)/B_ij - 4*DX(l,j,i,mu)*T_ij/pow(B_ij,2);
	}
	if (i!=mj) {
		number B_imj = a*a + DistanceSquared(l[i],l[mj]);
		res += -2.0*DX(l,i,pi,mu)/B_imj;
	}

	v[j*4+mu] += -f*res;
}

// mdV2r_nr
template <> void mdV2r_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& f, vec& v) {
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

	v[j*4+mu] += -f*res;
}

// ddV1r_nr
template <> void ddV1r_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
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
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// ddV2r_nr
template <> void ddV2r_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
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
	
	m(4*j+mu,4*k+nu) += f*res;
	
}

// mdI_nr<4>
template <> void mdI_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& f, vec& v) {
	if (mu==3) {
		uint nj = (j==0? (l.size()-1): j-1);
		v[j*4+mu] += -f*((l[nj])[2]-(l[j])[2]);
	}
	else if (mu==2) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		v[j*4+mu] += -f*((l[pj])[3]-(l[j])[3]);
	}
}

// ddI_nr<4>
template <> void ddI_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, \
							const Loop<4>& l, const number& f, mat& m) {
	if (mu==3 && nu==2) {
		uint nj = (j==0? (l.size()-1): j-1);
		if (k==nj)
			m(j*4+mu,k*4+nu) += f;
		else if (k==j)
			m(j*4+mu,k*4+nu) -= f;
	}
	else if (mu==2 && nu==3) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		if (k==pj)
			m(j*4+mu,k*4+nu) += f;
		else if (k==j)
			m(j*4+mu,k*4+nu) -= f;
	}
}	

// I0<4>
template <> void I0<4> (const uint& j, const Loop<4>& l, const number& f, number& result) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	result += f*(l[j])[2]*((l[pj])[3]-(l[j])[3]);
}					
		
