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

// FGamma
template <uint Dim>
void FGamma (const uint& j, const Loop<Dim>& l, const number& f, number& result) {

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

	for (uint i=0; i<l.size(); i++) {
		mdV1r_nr( j, mu, i, l, a, f, v)
	}
	
}

// mdV1r_nr
template<uint Dim>
void mdV1r_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v) {
	number denom, res = 0.0;
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	uint pi = (i==(l.size()-1)? 0: i+1);
	uint ni = (i==0? (l.size()-1): i-1);
		
	if (i!=j) {
		denom = a*a + DistanceSquared(l[i],l[j]);
		res -= 2.0*((l[pi])[mu]-(l[i])[mu])/denom;
		res -= 4.0*((l[j])[mu]-(l[i])[mu])*Dot(l[pi],l[i],l[pj],l[j])/pow(denom,2);
	}
	if (i!=nj) {
		denom = a*a + DistanceSquared(l[i],l[nj]);
		res += 2.0*((l[pi])[mu]-(l[i])[mu])/denom;
	}

	v[j*Dim+mu] += f*res;
}

// ddV1r_nr
template<uint Dim>
void ddV1r_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a const number& f, mat& m) {
	number res = 0.0;
	
	uint nj = (j==0? (l.size()-1): j-1);
	uint pj = (j==(l.size()-1)? 0: j+1);		
	uint nk = (k==0? (l.size()-1): k-1);
	uint pk = (k==(l.size()-1)? 0: k+1);
	
	number denom = a*a + DistanceSquared(l[j],l[k]);
	number denom_nj = a*a + DistanceSquared(l[nj],l[k]);
	number denom_nk = a*a + DistanceSquared(l[j],l[nk]);
	number denom_njk = a*a + DistanceSquared(l[nj],l[nk]);
	
	number orig_num = Dot(l[pj],l[j],l[pk],l[k]);
	number x_jk_mu = (l[j])[mu]-(l[k])[mu];
	number x_jk_nu = (l[j])[nu]-(l[k])[nu];
	number dx_j_mu = (l[pj])[mu]-(l[j])[mu];
	number dx_j_nu = (l[pj])[nu]-(l[j])[nu];
	number dx_k_mu = (l[pk])[mu]-(l[k])[mu];
	number dx_k_nu = (l[pk])[nu]-(l[k])[nu];
	
	// terms where mu==nu, without sums
	if (mu==nu) {
		if (k!=j)
			res += 2.0/denom_njk + 2.0/denom + 4.0*orig_num/pow(denom,2);
		if (k!=nj)
			res += -2.0/denom_nj;
		if (k!=pj)
			res += -2.0/denom_nk;
	}
	
	// terms where mu not nexcessarily equal to nu, without sums
	if (k!=j)
		res += 4.0*( -dx_k_mu*x_jk_nu + x_jk_mu*_dx_j_nu)/pow(denom,2)- 16.0*orig_num*x_jk_mu*x_jk_nu/pow(denom,3);
	if (k!=nj)
		res += 4.0*dx_k_mu*(l[nj])[nu]-(l[k])[nu])/pow(denom_nj,2);
	if (k!=pj)
		res += -4.0*dx_j_nu*(l[j])[mu]-(l[nk])[mu])/pow(denom_nk,2);
	
	// terms with sums
	if (k==j || k==nj || k==pj) {
	
		uint ni, pi;
		number dx_i_mu, dx_i_nu, denom_ij, denom_inj, orig_num_ij;
		
		for (uint i=0; i<l.size(); i++) {
		
			ni = (i==0? (l.size()-1): i-1);
			pi = (i==(l.size()-1)? 0: i+1);
			dx_i_mu = (l[pi])[mu]-(l[i])[mu];
			dx_i_nu = (l[pi])[nu]-(li])[nu];
			denom_ij = a*a + DistanceSquared(l[i],l[j]);
			denom_inj = a*a + DistanceSquared(l[i],l[nj]);
			
			if (k==j && i!=j)
				res += 	(4.0/pow(denom_inj,2)) * ( (l[j])[nu]-(l[i])[nu])*dx_i_mu + (l[j])[mu]-(l[i])[mu])*dx_i_nu \
						+ 2.0*((l[j])[mu]-(l[i])[mu])*((l[j])[nu]-(l[i])[nu])*Dot(l[pj],l[j],l[pi],l[i])/denom_ij );
			if (k==nj && i!=nj) 
				res += -4.0*(l[nj])[nu]-(l[i])[nu])*dx_i_mu/pow(denom_inj,2);
			if (k==pj && i!=pj) 
				res += -4.0*(l[j])[mu]-(l[i])[mu])*dx_i_nu/pow(denom_ij,2);
			
			if (mu==nu && i!=j)
				res += -4.0*Dot(l[pj],l[j],l[pi],l[i])/pow(denom_ij,2);
			
		}
		
	}
	
	m(Dim*j+mu,Dim*k+nu) += res;
	
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
template void FGamma<2>(const uint& j, const Loop<2>& l, const number& f, number& result);
template void mdL_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void ddL_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void mdS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& p, vec& v);
template void ddS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const number& p, mat& m);
template void mdsqrtS0_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& sqrt4s0, const number& p, vec& v);
template void ddsqrtS0_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l,\
								 const number& sqrt4s0, const number& p, mat& m);
template void mdV1r_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const number& a, const number& p, vec& v);
template void mdV1r_nr<2>(const uint& j, const uint& mu, const uint& i, const Loop<2>& l, const number& a, const number& p, vec& v);
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
template void FGamma<4>(const uint& j, const Loop<4>& l, const number& f, number& result);
template void mdL_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void ddL_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void mdS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& p, vec& v);
template void ddS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const number& p, mat& m);
template void mdsqrtS0_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& sqrt4s0, const number& p, vec& v);
template void ddsqrtS0_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l,\
								 const number& sqrt4s0, const number& p, mat& m);
template void mdV1r_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const number& a, const number& p, vec& v);
template void mdV1r_nr<4>(const uint& j, const uint& mu, const uint& i, const Loop<4>& l, const number& a, const number& p, vec& v);
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
		
