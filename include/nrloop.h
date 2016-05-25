/*
	header file for solving classical equations of motion for monopole worldline via newton-raphson method
*/
 
#ifndef __NRLOOP_H_INCLUDED__
#define __NRLOOP_H_INCLUDED__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include "simple.h"
#include "folder.h"
#include "parameters.h"
#include "genloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - typedefs
		1 - nr loop functions
		2 - loopToVector, vectorToLoop
		3 - filename functions
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0 - typedefs
----------------------------------------------------------------------------------------------------------------------------*/

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

/*----------------------------------------------------------------------------------------------------------------------------
	1 - nr loop functions
----------------------------------------------------------------------------------------------------------------------------*/

// L
template <uint Dim>
void L (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// LDisjoint
template <uint Dim>
void LDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// DistPow
template <uint Dim>
void DistPow (const uint& j, const Loop<Dim>& l, const number& w, const number& f, number& result);

// Sm
template <uint Dim>
void Sm (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// S0
template <uint Dim>
void S0 (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// S0Disjoint
template <uint Dim>
void S0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// I0
template <uint Dim>
void I0 (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// I0Disjoint
template <uint Dim>
void I0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// Angle
template <uint Dim>
void Angle (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// MaxXn
template <uint Dim>
void MaxXn (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result);

// PseudoAngle, in 3-4 plane
template <uint Dim>
void PseudoAngle (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// FGamma
template <uint Dim>
void FGamma (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// Vor
template <uint Dim>
void Vor (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Vlr
template <uint Dim>
void Vlr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Ver
template <uint Dim>
void Ver (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Vdr
template <uint Dim>
void Vdr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Vthr
template <uint Dim>
void Vthr (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// VthrDisjoint
template <uint Dim>
void VthrDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// Gaussian
template <uint Dim>
void Gaussian (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// GaussianDisjoint
template <uint Dim>
void GaussianDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// InlineCurvatureMax
template <uint Dim>
void InlineCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// InlineCurvatureMax
template <uint Dim>
void InlineCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// InlineCurvatureAvg
template <uint Dim>
void InlineCurvatureAvg (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// InlineCurvatureAvg
template <uint Dim>
void InlineCurvatureAvg (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// CuspCurvatureMax
template <uint Dim>
void CuspCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// KGAvg
template <uint Dim>
void KGAvg (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGAvg
template <uint Dim>
void KGAvg (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// KGAvgPlane
template <uint Dim>
void KGAvgPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGAvgPlane
template <uint Dim>
void KGAvgPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// mdPX_nr
template<uint Dim>
void mdPX_nr(const Loop<Dim>& l, const uint& loc, const Point<Dim>& P, const number& f, vec& v);

// mdL_nr
template<uint Dim>
void mdL_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v);

// ddL_nr
template<uint Dim>
void ddL_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdLDisjoint_nr
template<uint Dim>
void mdLDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);

// ddLDisjoint_nr
template<uint Dim>
void ddLDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
			const number& beta, const number& f, mat& m);

// mdDistPow_nr
template<uint Dim>
void mdDistPow_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& w, const number& f, vec& v);

// ddDistPow_nr
template<uint Dim>
void ddDistPow_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& w, const number& f, mat& m);

// mdI_nr
template<uint Dim>
void mdI_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v);

// ddI_nr
template<uint Dim>
void ddI_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdIDisjoint_nr
template<uint Dim>
void mdIDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);

// ddIDisjoint_nr
template<uint Dim>
void ddIDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& beta, const number& f, mat& m);

// mdS0_nr
template<uint Dim>
void mdS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v);

// ddS0_nr
template<uint Dim>
void ddS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdsqrtS0_nr
template<uint Dim>
void mdsqrtS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0, const number& f, vec& v);

// ddsqrtS0_nr
template<uint Dim>
void ddsqrtS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& f, mat& m);
						
// mdS0Disjoint_nr
template<uint Dim>
void mdS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);

// ddS0Disjoint_nr
template<uint Dim>
void ddS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdsqrtS0Disjoint_nr
template<uint Dim>
void mdsqrtS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0,\
 const number& beta, const number& f, vec& v);

// ddsqrtS0Disjoint_nr
template<uint Dim>
void ddsqrtS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& beta, const number& f, mat& m);
						
// mdVor_nr
template<uint Dim>
void mdVor_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVlr_nr
template<uint Dim>
void mdVlr_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVer_nr
template<uint Dim>
void mdVer_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVdr_nr
template<uint Dim>
void mdVdr_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVthr_nr
template<uint Dim>
void mdVthr_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v);

// mdVor_nr
template<uint Dim>
void mdVor_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVlr_nr
template<uint Dim>
void mdVlr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVer_nr
template<uint Dim>
void mdVer_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVdr_nr
template<uint Dim>
void mdVdr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVthr_nr
template<uint Dim>
void mdVthr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v);
					 	
// mdVthrDisjoint_nr
template<uint Dim>
void mdVthrDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v);

// mdGaussian_nr
template<uint Dim>
void mdGaussian_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdGaussianDisjoint_nr
template<uint Dim>
void mdGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);

// ddVor_nr
template<uint Dim>
void ddVor_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);

// ddVlr_nr
template<uint Dim>
void ddVlr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddVer_nr
template<uint Dim>
void ddVer_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddVdr_nr
template<uint Dim>
void ddVdr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddVthr_nr
template<uint Dim>
void ddVthr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddVthrDisjoint_nr
template<uint Dim>
void ddVthrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);

// ddGaussian_nr
template<uint Dim>
void ddGaussian_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddGaussianDisjoint_nr
template<uint Dim>
void ddGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// mdFGamma_nr
template<uint Dim>
void mdFGamma_nr(const Loop<Dim>& l, const uint& loc, const number& f, vec& v);

// ddFGamma_nr
template<uint Dim>
void ddFGamma_nr(const Loop<Dim>& l, const uint& loc, const number& f, mat& m);

/*----------------------------------------------------------------------------------------------------------------------------
	2 - loopToVector, vectorToLoop
----------------------------------------------------------------------------------------------------------------------------*/

// loopToVector
template<uint Dim>
void loopToVector(const Loop<Dim>& l, vec& v);

// vectorToLoop
template<uint Dim>
void vectorToLoop(const vec& v, Loop<Dim>& l);

// printAsLoop
void printAsLoop(const string& f, const uint& Dim, const vec& v);

// printAsLoop
void printAsLoop(const string& f, const uint& Dim, const vec& v, const uint len);

/*----------------------------------------------------------------------------------------------------------------------------
	3 - filename functions
----------------------------------------------------------------------------------------------------------------------------*/

// filenameLoop
template<uint Dim>
Filename filenameLoopNR(const Parameters&);

#endif // __NRLOOP_H_INCLUDED__
