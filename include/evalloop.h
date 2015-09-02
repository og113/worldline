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
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - loop functions
		- S0
		- V
----------------------------------------------------------------------------------------------------------------------------*/

// S0
template <uint Dim>
number S0 (const Loop<Dim>& l);

// V
template <uint Dim>
number V (const Loop<Dim>& l);


#endif // __EVALLOOP_H_INCLUDED__
