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

// In
template <uint Dim>
void In (const uint& j, const Loop<Dim>& l, const number& n, const number& f, number& result);

// InDisjoint
template <uint Dim>
void InDisjoint (const uint& j, const Loop<Dim>& l, const number& n, const number& beta, const number& f, number& result);

// Angle
template <uint Dim>
void Angle (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// MaxXn
template <uint Dim>
void MaxXn (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result);


// MinXnDisjoint
template <uint Dim>
void MinXnDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result);

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

// VthrDisjointLR
template <uint Dim>
void VthrDisjointLR (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);


// VnonrelDisjoint
template <uint Dim>
void VnonrelDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// VnonrelrDisjoint
template <uint Dim>
void VnonrelrDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// Gaussian
template <uint Dim>
void Gaussian (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// GaussianDisjoint
template <uint Dim>
void GaussianDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianLRDisjoint
template <uint Dim>
void GaussianLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianThermal
template <uint Dim>
void GaussianThermal (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// GaussianThermalDisjoint
template <uint Dim>
void GaussianThermalDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianThermalLRDisjoint
template <uint Dim>
void GaussianThermalLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianThermal2
template <uint Dim>
void GaussianThermal2 (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// GaussianThermal2Disjoint
template <uint Dim>
void GaussianThermal2Disjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianThermal2LRDisjoint
template <uint Dim>
void GaussianThermal2LRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// InlineCurvatureMax
template <uint Dim>
void InlineCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// InlineCurvatureMax
template <uint Dim>
void InlineCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// InlineCurvatureMaxDisjoint
template <uint Dim>
void InlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// InlineCurvatureMaxDisjoint
template <uint Dim>
void InlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);


// CuspCurvatureMax
template <uint Dim>
void CuspCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// CuspCurvatureMax
template <uint Dim>
void CuspCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// CuspCurvatureMaxDisjoint
template <uint Dim>
void CuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// CuspCurvatureMaxDisjoint
template <uint Dim>
void CuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGMax
template <uint Dim>
void KGMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGMaxPlane
template <uint Dim>
void KGMaxPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// KGMaxDisjoint
template <uint Dim>
void KGMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// KGMaxDisjoint
template <uint Dim>
void KGMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);

// KGMaxPlaneDisjoint
template <uint Dim>
void KGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// KGMaxPlaneDisjoint
template <uint Dim>
void KGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);

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

// mdI0_nr
template<uint Dim>
void mdI0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v);

// ddI0_nr
template<uint Dim>
void ddI0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdI0Disjoint_nr
template<uint Dim>
void mdI0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);

// ddI0Disjoint_nr
template<uint Dim>
void ddI0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& beta, const number& f, mat& m);

// mdIn_nr
template<uint Dim>
void mdIn_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& f, vec& v);

// ddIn_nr
template<uint Dim>
void ddIn_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& f, mat& m);

// mdInDisjoint_nr
template<uint Dim>
void mdInDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, vec& v);

// ddInDisjoint_nr
template<uint Dim>
void ddInDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, mat& m);

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
	
// mdVthrDisjointLR_nr
template<uint Dim>
void mdVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v);
				 	
// mdVnonrelDisjoint_nr
template<uint Dim>
void mdVnonrelDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);
					 	
// mdVnonrelrDisjoint_nr
template<uint Dim>
void mdVnonrelrDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& a\
						, const number& f, vec& v);

// mdGaussian_nr
template<uint Dim>
void mdGaussian_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdGaussianDisjoint_nr
template<uint Dim>
void mdGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianLRDisjoint_nr
template<uint Dim>
void mdGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermal_nr
template<uint Dim>
void mdGaussianThermal_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v);

// mdGaussianThermalDisjoint_nr
template<uint Dim>
void mdGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermalLRDisjoint_nr
template<uint Dim>
void mdGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermal2_nr
template<uint Dim>
void mdGaussianThermal2_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v);

// mdGaussianThermal2Disjoint_nr
template<uint Dim>
void mdGaussianThermal2Disjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermal2LRDisjoint_nr
template<uint Dim>
void mdGaussianThermal2LRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
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

// ddVthrDisjointLR_nr
template<uint Dim>
void ddVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);

// ddVnonrelDisjoint_nr
template<uint Dim>
void ddVnonrelDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& f, mat& m);
						 
// ddVnonrelrDisjoint_nr
template<uint Dim>
void ddVnonrelrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);

// ddGaussian_nr
template<uint Dim>
void ddGaussian_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddGaussianDisjoint_nr
template<uint Dim>
void ddGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m); 

// ddGaussianLRDisjoint_nr
template<uint Dim>
void ddGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m); 

// ddGaussianThermal_nr
template<uint Dim>
void ddGaussianThermal_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermalDisjoint_nr
template<uint Dim>
void ddGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermalLRDisjoint_nr
template<uint Dim>
void ddGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
	
// ddGaussianThermal2_nr
template<uint Dim>
void ddGaussianThermal2_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermal2Disjoint_nr
template<uint Dim>
void ddGaussianThermal2Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermal2LRDisjoint_nr
template<uint Dim>
void ddGaussianThermal2LRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
					 
// mdFGamma_nr
template<uint Dim>
void mdFGamma_nr(const Loop<Dim>& l, const uint& loc, const number& f, vec& v);

// ddFGamma_nr
template<uint Dim>
void ddFGamma_nr(const Loop<Dim>& l, const uint& loc, const number& f, mat& m);

// ErgS0_nr
template<uint Dim>
void ErgS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, number& erg);

// ErgS0Disjoint_nr
template<uint Dim>
void ErgS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);

// ErgL_nr
template<uint Dim>
void ErgL_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, number& erg);

// ErgLDisjoint_nr
template<uint Dim>
void ErgLDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);

// ErgsqrtS0_nr
template<uint Dim>
void ErgsqrtS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, number& erg);

// ErgsqrtS0Disjoint_nr
template<uint Dim>
void ErgsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, number& erg);

// ErgI0_nr
template<uint Dim>
void ErgI0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, number& erg);

// ErgI0Disjoint_nr
template<uint Dim>
void ErgI0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);

// ErgIn_nr
template<uint Dim>
void ErgIn_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& f, number& erg);

// ErgInDisjoint_nr
template<uint Dim>
void ErgInDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, number& erg);

// ErgVor_nr
template<uint Dim>
void ErgVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg);

// ErgVthr_nr
template<uint Dim>
void ErgVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgVthrDisjoint_nr
template<uint Dim>
void ErgVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgVthrDisjointLR_nr
template<uint Dim>
void ErgVthrDisjointLR_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussian_nr
template<uint Dim>
void ErgGaussian_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, number& erg);

// ErgGaussianDisjoint_nr
template<uint Dim>
void ErgGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianLRDisjoint_nr
template<uint Dim>
void ErgGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermal_nr
template<uint Dim>
void ErgGaussianThermal_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermalDisjoint_nr
template<uint Dim>
void ErgGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermalLRDisjoint_nr
template<uint Dim>
void ErgGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermal2_nr
template<uint Dim>
void ErgGaussianThermal2_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermal2Disjoint_nr
template<uint Dim>
void ErgGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianTherma2lLRDisjoint_nr
template<uint Dim>
void ErgGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// PS0_nr
template<uint Dim>
void PS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v);

// PS0Disjoint_nr
template<uint Dim>
void PS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);

// PL_nr
template<uint Dim>
void PL_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v);

// PLDisjoint_nr
template<uint Dim>
void PLDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);

// PsqrtS0_nr
template<uint Dim>
void PsqrtS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, vec& v);

// PsqrtS0Disjoint_nr
template<uint Dim>
void PsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, vec& v);

// PI0_nr
template<uint Dim>
void PI0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v);

// PI0Disjoint_nr
template<uint Dim>
void PI0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);

// PIn_nr
template<uint Dim>
void PIn_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& f, vec& v);

// PInDisjoint_nr
template<uint Dim>
void PInDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, vec& v);

// PVor_nr
template<uint Dim>
void PVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v);

// PRVor_nr
template<uint Dim>
void PRVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);

// PVthr_nr
template<uint Dim>
void PVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRVthr_nr
template<uint Dim>
void PRVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PVthrDisjoint_nr
template<uint Dim>
void PVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRVthrDisjoint_nr
template<uint Dim>
void PRVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussian_nr
template<uint Dim>
void PGaussian_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, vec& v);

// PRGaussian_nr
template<uint Dim>
void PRGaussian_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);

// PGaussianDisjoint_nr
template<uint Dim>
void PGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianDisjoint_nr
template<uint Dim>
void PRGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianLRDisjoint_nr
template<uint Dim>
void PGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianLRDisjoint_nr
template<uint Dim>
void PRGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermal_nr
template<uint Dim>
void PGaussianThermal_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal_nr
template<uint Dim>
void PRGaussianThermal_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermalDisjoint_nr
template<uint Dim>
void PGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermalDisjoint_nr
template<uint Dim>
void PRGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermalLRDisjoint_nr
template<uint Dim>
void PGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermalLRDisjoint_nr
template<uint Dim>
void PRGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal2_nr
template<uint Dim>
void PRGaussianThermal2_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermal2Disjoint_nr
template<uint Dim>
void PGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal2Disjoint_nr
template<uint Dim>
void PRGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermal2LRDisjoint_nr
template<uint Dim>
void PGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal2LRDisjoint_nr
template<uint Dim>
void PRGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

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

// dimReduce
void dimReduce(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout);

// dimIncrease
void dimIncrease(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout);

/*----------------------------------------------------------------------------------------------------------------------------
	3 - filename functions
----------------------------------------------------------------------------------------------------------------------------*/

// filenameLoop
template<uint Dim>
Filename filenameLoopNR(const Parameters&, const string&);

// filenameLoop
template<uint Dim>
Filename filenameThermalNR(const Parameters&, const string&);

#endif // __NRLOOP_H_INCLUDED__
