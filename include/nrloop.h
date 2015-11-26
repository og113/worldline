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
#include "parameters.h"
#include "genloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - typedefs
		1 - nr vector and matrix generation functions
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0 - typedefs
----------------------------------------------------------------------------------------------------------------------------*/

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

/*----------------------------------------------------------------------------------------------------------------------------
	1 - nr problem generation functions
----------------------------------------------------------------------------------------------------------------------------*/

// dL_nr
template<uint Dim>
void dL_nr(const uint& j, const uint& mu, const Loop<Dim>& l, vec& v);

// ddL_nr
template<uint Dim>
void ddL_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, mat& m);

// dI_nr
template<uint Dim>
void dI_nr(const uint& j, const uint& mu, const Loop<Dim>& l, vec& v);

// ddI_nr
template<uint Dim>
void ddI_nr(const uint& j, const uint& mu, const uint& k, const uint& nu, const Loop<Dim>& l, mat& m);

#endif // __NRLOOP_H_INCLUDED__
