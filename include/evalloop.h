/*
	header file for calculating quantities from unit loops in D dimensions for use in worldline programs.
*/
 
#ifndef __EVALLOOP_H_INCLUDED__
#define __EVALLOOP_H_INCLUDED__

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_exp.h> // exponential function
#include "genloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - loop functions
		2 - data struct
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - loop functions
		- S0
		- DS0
		- V0
		- aprxDV0
----------------------------------------------------------------------------------------------------------------------------*/

// S0
template <uint Dim>
number S0 (const Loop<Dim>& l);

// DS0
template <uint Dim>
number DS0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

// V0
template <uint Dim>
number V0 (const Loop<Dim>& l);

// aprxDV0
template <uint Dim>
number aprxDV0 (const Loop<Dim>& l, const Point<Dim>& p, const uint& loc);

/*----------------------------------------------------------------------------------------------------------------------------
	2 - data struct
----------------------------------------------------------------------------------------------------------------------------*/

// data


#endif // __EVALLOOP_H_INCLUDED__
