/*
	header file for solving classical equations of motion for monopole worldline via newton-raphson method
*/
 
#ifndef __NRLOOP2_H_INCLUDED__
#define __NRLOOP2_H_INCLUDED__

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
		2 - NEWloopToVector, vectorToLoop
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
void NEWL (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// LDisjoint
template <uint Dim>
void NEWLDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// DistPow
template <uint Dim>
void NEWDistPow (const uint& j, const Loop<Dim>& l, const number& w, const number& f, number& result);

// Sm
template <uint Dim>
void NEWSm (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// S0
template <uint Dim>
void NEWS0 (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// S0Disjoint
template <uint Dim>
void NEWS0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// I0
template <uint Dim>
void NEWI0 (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// I0Disjoint
template <uint Dim>
void NEWI0Disjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// In
template <uint Dim>
void NEWIn (const uint& j, const Loop<Dim>& l, const number& n, const number& f, number& result);

// InDisjoint
template <uint Dim>
void NEWInDisjoint (const uint& j, const Loop<Dim>& l, const number& n, const number& beta, const number& f, number& result);

// Angle
template <uint Dim>
void NEWAngle (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// MaxXn
template <uint Dim>
void NEWMaxXn (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result);


// MinXnDisjoint
template <uint Dim>
void NEWMinXnDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const uint& n, const number& f, number& result);

// NEWPseudoAngle, in 3-4 plane
template <uint Dim>
void NEWPseudoAngle (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// FGamma
template <uint Dim>
void NEWFGamma (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// Vor
template <uint Dim>
void NEWVor (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Vlr
template <uint Dim>
void NEWVlr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Ver
template <uint Dim>
void NEWVer (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Vdr
template <uint Dim>
void NEWVdr (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// Vthr
template <uint Dim>
void NEWVthr (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// VthrDisjoint
template <uint Dim>
void NEWVthrDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// VthrDisjointLR
template <uint Dim>
void NEWVthrDisjointLR (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);


// VnonrelDisjoint
template <uint Dim>
void NEWVnonrelDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// VnonrelrDisjoint
template <uint Dim>
void NEWVnonrelrDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// Gaussian
template <uint Dim>
void NEWGaussian (const uint& j, const uint& k, const Loop<Dim>& l, const number& a, const number& f, number& result);

// GaussianDisjoint
template <uint Dim>
void NEWGaussianDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianLRDisjoint
template <uint Dim>
void NEWGaussianLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianThermal
template <uint Dim>
void NEWGaussianThermal (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// GaussianThermalDisjoint
template <uint Dim>
void NEWGaussianThermalDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// GaussianThermalLRDisjoint
template <uint Dim>
void NEWGaussianThermalLRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// NEWGaussianThermal2
template <uint Dim>
void NEWGaussianThermal2 (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a, const number& f, number& result);

// NEWGaussianThermal2Disjoint
template <uint Dim>
void NEWGaussianThermal2Disjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// NEWGaussianThermal2LRDisjoint
template <uint Dim>
void NEWGaussianThermal2LRDisjoint (const uint& j, const uint& k, const Loop<Dim>& l, const number& beta, const number& a,const number& f, number& result);

// InlineCurvatureMax
template <uint Dim>
void NEWInlineCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// InlineCurvatureMax
template <uint Dim>
void NEWInlineCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// InlineCurvatureMaxDisjoint
template <uint Dim>
void NEWInlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// InlineCurvatureMaxDisjoint
template <uint Dim>
void NEWInlineCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);


// CuspCurvatureMax
template <uint Dim>
void NEWCuspCurvatureMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// CuspCurvatureMax
template <uint Dim>
void NEWCuspCurvatureMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// CuspCurvatureMaxDisjoint
template <uint Dim>
void NEWCuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// CuspCurvatureMaxDisjoint
template <uint Dim>
void NEWCuspCurvatureMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);

// KGMax
template <uint Dim>
void NEWKGMax (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGMax
template <uint Dim>
void NEWKGMax (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// KGMaxPlane
template <uint Dim>
void NEWKGMaxPlane (const uint& j, const Loop<Dim>& l, const number& f, number& result);

// KGMaxPlane
template <uint Dim>
void NEWKGMaxPlane (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& f, number& result);

// KGMaxDisjoint
template <uint Dim>
void NEWKGMaxDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// KGMaxDisjoint
template <uint Dim>
void NEWKGMaxDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);

// KGMaxPlaneDisjoint
template <uint Dim>
void NEWKGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const number& beta, const number& f, number& result);

// KGMaxPlaneDisjoint
template <uint Dim>
void NEWKGMaxPlaneDisjoint (const uint& j, const Loop<Dim>& l, const uint& ex1, const uint& ex2, const number& beta, const number& f, number& result);

// mdPX_nr
template<uint Dim>
void NEWmdPX_nr(const Loop<Dim>& l, const uint& loc, const Point<Dim>& P, const number& f, vec& v);

// mdL_nr
template<uint Dim>
void NEWmdL_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v);

// ddL_nr
template<uint Dim>
void NEWddL_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdLDisjoint_nr
template<uint Dim>
void NEWmdLDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);

// ddLDisjoint_nr
template<uint Dim>
void NEWddLDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
			const number& beta, const number& f, mat& m);

// mdDistPow_nr
template<uint Dim>
void NEWmdDistPow_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& w, const number& f, vec& v);

// ddDistPow_nr
template<uint Dim>
void NEWddDistPow_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& w, const number& f, mat& m);

// mdI0_nr
template<uint Dim>
void NEWmdI0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v);

// ddI0_nr
template<uint Dim>
void NEWddI0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdI0Disjoint_nr
template<uint Dim>
void NEWmdI0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);

// ddI0Disjoint_nr
template<uint Dim>
void NEWddI0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& beta, const number& f, mat& m);

// mdIn_nr
template<uint Dim>
void NEWmdIn_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& f, vec& v);

// ddIn_nr
template<uint Dim>
void NEWddIn_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& f, mat& m);

// mdInDisjoint_nr
template<uint Dim>
void NEWmdInDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, vec& v);

// ddInDisjoint_nr
template<uint Dim>
void NEWddInDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& n, const number& beta, const number& f, mat& m);

// mdS0_nr
template<uint Dim>
void NEWmdS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& f, vec& v);

// ddS0_nr
template<uint Dim>
void NEWddS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdsqrtS0_nr
template<uint Dim>
void NEWmdsqrtS0_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0, const number& f, vec& v);

// ddsqrtS0_nr
template<uint Dim>
void NEWddsqrtS0_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& f, mat& m);
						
// mdS0Disjoint_nr
template<uint Dim>
void NEWmdS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);

// ddS0Disjoint_nr
template<uint Dim>
void NEWddS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, const number& f, mat& m);

// mdsqrtS0Disjoint_nr
template<uint Dim>
void NEWmdsqrtS0Disjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& sqrt4s0,\
 const number& beta, const number& f, vec& v);

// ddsqrtS0Disjoint_nr
template<uint Dim>
void NEWddsqrtS0Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, \
						const Loop<Dim>& l, const number& sqrt4s0, const number& beta, const number& f, mat& m);

// mdVor_nr
template<uint Dim>
void NEWmdVor_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVlr_nr
template<uint Dim>
void NEWmdVlr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVer_nr
template<uint Dim>
void NEWmdVer_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVdr_nr
template<uint Dim>
void NEWmdVdr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdVthr_nr
template<uint Dim>
void NEWmdVthr_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v);
					 	
// mdVthrDisjoint_nr
template<uint Dim>
void NEWmdVthrDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v);
	
// mdVthrDisjointLR_nr
template<uint Dim>
void NEWmdVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
					 	const number& beta, const number& a, const number& f, vec& v);
				 	
// mdVnonrelDisjoint_nr
template<uint Dim>
void NEWmdVnonrelDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& f, vec& v);
					 	
// mdVnonrelrDisjoint_nr
template<uint Dim>
void NEWmdVnonrelrDisjoint_nr(const uint& j, const uint& mu, const Loop<Dim>& l, const number& beta, const number& a\
						, const number& f, vec& v);

// mdGaussian_nr
template<uint Dim>
void NEWmdGaussian_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& a, const number& f, vec& v);

// mdGaussianDisjoint_nr
template<uint Dim>
void NEWmdGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianLRDisjoint_nr
template<uint Dim>
void NEWmdGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermal_nr
template<uint Dim>
void NEWmdGaussianThermal_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v);

// mdGaussianThermalDisjoint_nr
template<uint Dim>
void NEWmdGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWmdGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermal2_nr
template<uint Dim>
void NEWmdGaussianThermal2_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l, const number& beta, const number& a, const number& f, vec& v);

// mdGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWmdGaussianThermal2Disjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
			
// mdGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWmdGaussianThermal2LRDisjoint_nr(const uint& j, const uint& mu, const uint& i, const Loop<Dim>& l,\
			const number& beta, const number& a, const number& f, vec& v);
// ddVor_nr
template<uint Dim>
void NEWddVor_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);

// ddVlr_nr
template<uint Dim>
void NEWddVlr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddVer_nr
template<uint Dim>
void NEWddVer_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddVdr_nr
template<uint Dim>
void NEWddVdr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddVthr_nr
template<uint Dim>
void NEWddVthr_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddVthrDisjoint_nr
template<uint Dim>
void NEWddVthrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);

// ddVthrDisjointLR_nr
template<uint Dim>
void NEWddVthrDisjointLR_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);

// ddVnonrelDisjoint_nr
template<uint Dim>
void NEWddVnonrelDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& f, mat& m);
						 
// ddVnonrelrDisjoint_nr
template<uint Dim>
void NEWddVnonrelrDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);

// ddGaussian_nr
template<uint Dim>
void NEWddGaussian_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& a, const number& f, mat& m);
						 
// ddGaussianDisjoint_nr
template<uint Dim>
void NEWddGaussianDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m); 

// ddGaussianLRDisjoint_nr
template<uint Dim>
void NEWddGaussianLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m); 

// ddGaussianThermal_nr
template<uint Dim>
void NEWddGaussianThermal_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermalDisjoint_nr
template<uint Dim>
void NEWddGaussianThermalDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWddGaussianThermalLRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
	
// ddGaussianThermal2_nr
template<uint Dim>
void NEWddGaussianThermal2_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWddGaussianThermal2Disjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
						 
// ddGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWddGaussianThermal2LRDisjoint_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l,\
						 const number& beta, const number& a, const number& f, mat& m);
					 
// mdFGamma_nr
template<uint Dim>
void NEWmdFGamma_nr(const Loop<Dim>& l, const uint& loc, const number& f, vec& v);

// ddFGamma_nr
template<uint Dim>
void NEWddFGamma_nr(const Loop<Dim>& l, const uint& loc, const number& f, mat& m);

// ErgS0_nr
template<uint Dim>
void NEWErgS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, number& erg);

// ErgS0Disjoint_nr
template<uint Dim>
void NEWErgS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);

// ErgL_nr
template<uint Dim>
void NEWErgL_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, number& erg);

// ErgLDisjoint_nr
template<uint Dim>
void NEWErgLDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);

// ErgsqrtS0_nr
template<uint Dim>
void NEWErgsqrtS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, number& erg);

// ErgsqrtS0Disjoint_nr
template<uint Dim>
void NEWErgsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, number& erg);

// ErgI0_nr
template<uint Dim>
void NEWErgI0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, number& erg);

// ErgI0Disjoint_nr
template<uint Dim>
void NEWErgI0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, number& erg);

// ErgIn_nr
template<uint Dim>
void NEWErgIn_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& f, number& erg);

// ErgInDisjoint_nr
template<uint Dim>
void NEWErgInDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, number& erg);

// ErgVor_nr
template<uint Dim>
void NEWErgVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, number& erg);

// ErgVthr_nr
template<uint Dim>
void NEWErgVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgVthrDisjoint_nr
template<uint Dim>
void NEWErgVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgVthrDisjointLR_nr
template<uint Dim>
void NEWErgVthrDisjointLR_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussian_nr
template<uint Dim>
void NEWErgGaussian_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, number& erg);

// ErgGaussianDisjoint_nr
template<uint Dim>
void NEWErgGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianLRDisjoint_nr
template<uint Dim>
void NEWErgGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermal_nr
template<uint Dim>
void NEWErgGaussianThermal_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermalDisjoint_nr
template<uint Dim>
void NEWErgGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWErgGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermal2_nr
template<uint Dim>
void NEWErgGaussianThermal2_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWErgGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// ErgGaussianTherma2lLRDisjoint_nr
template<uint Dim>
void NEWErgGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, number& erg);

// PS0_nr
template<uint Dim>
void NEWPS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v);

// PS0Disjoint_nr
template<uint Dim>
void NEWPS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);

// PL_nr
template<uint Dim>
void NEWPL_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v);

// PLDisjoint_nr
template<uint Dim>
void NEWPLDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);

// PsqrtS0_nr
template<uint Dim>
void NEWPsqrtS0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& f, vec& v);

// PsqrtS0Disjoint_nr
template<uint Dim>
void NEWPsqrtS0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& sqrt4s0, const number& beta, const number& f, vec& v);

// PI0_nr
template<uint Dim>
void NEWPI0_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& f, vec& v);

// PI0Disjoint_nr
template<uint Dim>
void NEWPI0Disjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& beta, const number& f, vec& v);

// PIn_nr
template<uint Dim>
void NEWPIn_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& f, vec& v);

// PInDisjoint_nr
template<uint Dim>
void NEWPInDisjoint_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const number& n, const number& beta, const number& f, vec& v);

// PVor_nr
template<uint Dim>
void NEWPVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& a, const number& f, vec& v);

// PRVor_nr
template<uint Dim>
void NEWPRVor_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);

// PVthr_nr
template<uint Dim>
void NEWPVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRVthr_nr
template<uint Dim>
void NEWPRVthr_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PVthrDisjoint_nr
template<uint Dim>
void NEWPVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRVthrDisjoint_nr
template<uint Dim>
void NEWPRVthrDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussian_nr
template<uint Dim>
void NEWPGaussian_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& a, const number& f, vec& v);

// PRGaussian_nr
template<uint Dim>
void NEWPRGaussian_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& a, const number& f, vec& v);

// PGaussianDisjoint_nr
template<uint Dim>
void NEWPGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianDisjoint_nr
template<uint Dim>
void NEWPRGaussianDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianLRDisjoint_nr
template<uint Dim>
void NEWPGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianLRDisjoint_nr
template<uint Dim>
void NEWPRGaussianLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermal_nr
template<uint Dim>
void NEWPGaussianThermal_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal_nr
template<uint Dim>
void NEWPRGaussianThermal_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermalDisjoint_nr
template<uint Dim>
void NEWPGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermalDisjoint_nr
template<uint Dim>
void NEWPRGaussianThermalDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWPGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermalLRDisjoint_nr
template<uint Dim>
void NEWPRGaussianThermalLRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal2_nr
template<uint Dim>
void NEWPRGaussianThermal2_nr(const Loop<Dim>& l, const uint& loc, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWPGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal2Disjoint_nr
template<uint Dim>
void NEWPRGaussianThermal2Disjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

// PGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWPGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const number& beta, const number& a, const number& f, vec& v);

// PRGaussianThermal2LRDisjoint_nr
template<uint Dim>
void NEWPRGaussianThermal2LRDisjoint_nr(const Loop<Dim>& l, const uint& j, const uint& mu, const uint& k, const uint& nu, const number& beta, const number& a, const number& f, vec& v);

/*----------------------------------------------------------------------------------------------------------------------------
	2 - NEWloopToVector, vectorToLoop
----------------------------------------------------------------------------------------------------------------------------*/

// loopToVector
template<uint Dim>
void NEWloopToVector(const Loop<Dim>& l, vec& v);

// vectorToLoop
template<uint Dim>
void NEWvectorToLoop(const vec& v, Loop<Dim>& l);

// printAsLoop
void NEWprintAsLoop(const string& f, const uint& Dim, const vec& v);

// printAsLoop
void NEWprintAsLoop(const string& f, const uint& Dim, const vec& v, const uint len);

// dimReduce
void NEWdimReduce(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout);

// dimIncrease
void NEWdimIncrease(const vec& vin, const uint& dimin, const uint& Nin, vec& vout, const uint& dimout, const uint& zmout);

/*----------------------------------------------------------------------------------------------------------------------------
	3 - filename functions
----------------------------------------------------------------------------------------------------------------------------*/

// filenameLoop
template<uint Dim>
Filename NEWfilenameLoopNR(const Parameters&, const string&);

// filenameLoop
template<uint Dim>
Filename NEWfilenameThermalNR(const Parameters&, const string&);

#endif // __NRLOOP2_H_INCLUDED__
