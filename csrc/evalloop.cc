/*
	definition file for calculating quantities from unit loops in D dimensions for use in worldline programs.
*/

#include "evalloop.h"


/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - loop functions
		2 - explicit instantiation
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - loop functions
		- S0
		- V0
		
n.b. functions are defined for unit loops. for other loops the result must be multiplied by appropriate factors (depending on Dim) of T, the length of the loop.
----------------------------------------------------------------------------------------------------------------------------*/

// S0
template <uint Dim>
number S0 (const Loop<Dim>& l) {
	number result = DistanceSquared(l[l.size()-1],l[0]);
	for (uint j=0; j<(l.size()-1); j++)
		result += DistanceSquared(l[j+1],l[j]);
	return result*l.size()/4.0;
}

// V0
template <uint Dim>
number V0 (const Loop<Dim>& l) {
	number result = 0.0;
	for (uint j=0; j<l.size(); j++) {
		for (uint k=0; k<l.size(); k++) {
			if (j!=k)
				result += pow(DistanceSquared(l[j],l[k]),(2.0-Dim)/2.0);
		}
	}
	return result/pow(l.size(),2);
}


/*----------------------------------------------------------------------------------------------------------------------------
	4 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------*/
// Dim=4
template number S0<4> (const Loop<4>& l);

// V0, Dim=4, slightly changed for speed
template <> number V0 <4>(const Loop<4>& l) {
	number result = 0.0;
	for (uint j=0; j<l.size(); j++) {
		for (uint k=0; k<l.size(); k++) {
			if (j!=k)
				result += 1.0/DistanceSquared(l[j],l[k]);
		}
	}
	return result/pow(l.size(),2);
}

// Dim=2
template number S0<2> (const Loop<2>& l);
template number V0<2> (const Loop<2>& l);
