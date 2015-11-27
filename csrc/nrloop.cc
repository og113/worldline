/*
	definitions for solving classical equations of motion for monopole worldline via newton-raphson method
*/

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include "simple.h"
#include "parameters.h"
#include "genloop.h"
#include "nrloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - nr vector and matrix generation functions
		2 - explicit instatiation
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - nr problem generation functions
----------------------------------------------------------------------------------------------------------------------------*/

// dL_nr
template<uint Dim>
void dL_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const Parameters& p, vec& v) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	number norm = Distance(l[j],l[nj]), temp;
	temp = ((l[j])[mu]-(l[nj])[mu])/norm;
	norm = Distance(l[pj],l[j]);
	temp += ((l[pj])[mu]-(l[j])[mu])/norm;
	v[j*Dim+mu] += temp;
}

// ddL_nr
template<uint Dim>
void ddL_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const Parameters& p, mat& m) {
	if (k==j) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		uint nj = (j==0? (l.size()-1): j-1);
		number normn = Distance(l[j],l[nj]), normp = Distance(l[pj],l[j]), temp;
		if (mu==nu)
			m(j*Dim+mu,j*Dim+mu) += 1.0/normn + 1.0/normp;
		temp = -((l[j])[mu]-(l[nj])[mu])*((l[j])[nu]-(l[nj])[nu])/pow(normn,3);
		temp -=((l[j])[mu]-(l[pj])[mu])*((l[j])[nu]-(l[pj])[nu])/pow(normp,3);
		m(j*Dim+mu,j*Dim+nu) += temp;
	}
	else if (k==(j-1)) {
		uint nj = (j==0? (l.size()-1): j-1);
		number norm = Distance(l[j],l[nj]), temp;
		if (mu==nu)
			m(j*Dim+mu,k*Dim+mu) -= 1.0/norm;
		temp = ((l[j])[mu]-(l[nj])[mu])*((l[j])[nu]-(l[nj])[nu])/pow(norm,3);
		m(j*Dim+mu,k*Dim+nu) += temp;
	}
	else if (k==(j+1)) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		number norm = Distance(l[j],l[pj]), temp;
		if (mu==nu)
			m(j*Dim+mu,k*Dim+mu) -= 1.0/norm;
		temp = ((l[j])[mu]-(l[pj])[mu])*((l[j])[nu]-(l[pj])[nu])/pow(norm,3);
		m(j*Dim+mu,k*Dim+nu) += temp;
	}
}

// dI_nr
template<uint Dim>
void dI_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const Parameters& p, vec& v) {
	if (mu==3) {
		uint nj = (j==0? (l.size()-1): j-1);
		v[j*Dim+mu] += (l[nj])[2]-(l[j])[2];
	}
	else if (mu==2) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		v[j*Dim+mu] += (l[pj])[3]-(l[j])[3];
	}
}

// ddI_nr
template<uint Dim>
void ddI_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const Parameters& p, mat& m) {
	if (mu==3 && nu==2) {
		uint pj = (j==(l.size()-1)? 0: j+1);
		if (k==pj)
			m(j*Dim+mu,k*Dim+nu) += 1.0;
		else if (k==j)
			m(j*Dim+mu,k*Dim+nu) -= 1.0;
	}
	else if (mu==2 && nu==3) {
		uint nj = (j==0? (l.size()-1): j-1);
		if (k==nj)
			m(j*Dim+mu,k*Dim+nu) += 1.0;
		else if (k==j)
			m(j*Dim+mu,k*Dim+nu) -= 1.0;
	}
}

// dS0_nr
template<uint Dim>
void dS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const Parameters& p, vec& v) {
	uint pj = (j==(l.size()-1)? 0: j+1);
	uint nj = (j==0? (l.size()-1): j-1);
	v[j*Dim+mu] += (2.0*(l[j])[mu] - (l[nj])[mu] - (l[pj])[mu])/2.0;
}

// ddS0_nr
template<uint Dim>
void ddS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, mat& m) {
	if (mu==nu) {
		if (k==j) {
			m(j*Dim+mu,j*Dim+mu) += 1.0;
		}
		else if (k==(j-1)) {
			m(j*Dim+mu,k*Dim+mu) -= 0.5;
		}
		else if (k==(j+1)) {
			m(j*Dim+mu,k*Dim+mu) -= 0.5;
		}
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	2 - explicit instantiation
----------------------------------------------------------------------------------------------------------------------------*/

// dim 2
template void dI_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const Parameters& p, vec& v);
template void ddI_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const Parameters& p, mat& m);
template void dL_nr<2>(const uint& j, const uint& mu, const Loop<2>& l, const Parameters& p, vec& v);
template void ddL_nr<2>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<2>& l, const Parameters& p, mat& m);

// dim 4
template void dI_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const Parameters& p, vec& v);
template void ddI_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const Parameters& p, mat& m);
template void dL_nr<4>(const uint& j, const uint& mu, const Loop<4>& l, const Parameters& p, vec& v);
template void ddL_nr<4>(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<4>& l, const Parameters& p, mat& m);
