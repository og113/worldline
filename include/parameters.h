/*-------------------------------------------------------------------------------------------------------------------------
	parameters for worldline calculations
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __PARAMETERS_H_INCLUDED__
#define __PARAMETERS_H_INCLUDED__

#include <iostream>
#include <vector>
#include "simple.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - Parameters
	2 - ParametersRange
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

// ParametersRange declaration
struct ParametersRange;

// Parameters
struct Parameters {
	static const uint Size;
	enum Label { nl=1, ng=2, nms=3, k=4, g=5};
	uint Nl;
	uint Ng;
	uint Nms;
	uint K;
	number G;
	void step(const ParametersRange&);
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


/*-------------------------------------------------------------------------------------------------------------------------
	2. ParametersRange
		- ParametersRange
		- operator<<
		- operator==
	
	N.B. for structs the compiler automatically writes the following member functions, unless user provided:
	 empty constructor; copy constructor; assignment operator=; destructor
-------------------------------------------------------------------------------------------------------------------------*/

// ParametersRange
struct ParametersRange {
	static const uint Size;
	ParametersRange();
	ParametersRange(const Parameters& min, const Parameters& max, const vector<uint>& steps);
	Parameters 		Min;
	Parameters 		Max;
	vector<uint> 	Steps;
	bool			toStep(Parameters::Label&) const;
	void save(const string& filename) const;
	void load(const string& filename);
	bool empty() const;
	ostream& writeBinary(ostream&) const;
	istream& readBinary(istream&);
};

// operator<<
ostream& operator<<(ostream&, const ParametersRange&);

// operator==
bool operator==(const Parameters& lhs, const ParametersRange& rhs);

#endif // __PARAMETERS_H_INCLUDED__
