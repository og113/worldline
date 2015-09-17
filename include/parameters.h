/*-------------------------------------------------------------------------------------------------------------------------
	parameters for worldline calculations
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __PARAMETERS_H_INCLUDED__
#define __PARAMETERS_H_INCLUDED__

#include <iostream>
#include "simple.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - Parameters
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. Parameters
		- Parameters
		- operator<<
		- operator==
	
	N.B. for structs the compiler automatically writes the following member functions, unless user provided:
	 empty constructor; copy constructor; assignment operator=; destructor
-------------------------------------------------------------------------------------------------------------------------*/

// Parameters
struct Parameters {
	uint LoopMin;
	uint LoopMax;
	uint Ng;
	uint Nms;
	uint K;
	number g;
	void save(const string& filename) const;
	void load(const string& filename);
	bool empty() const;
	ostream& writeBinary(ostream&) const;
	istream& readBinary(istream&);
};

// operator<<
ostream& operator<<(ostream&, const Parameters&);

// operator==
bool operator==(const Parameters& lhs, const Parameters& rhs);

#endif // __PARAMETERS_H_INCLUDED__
