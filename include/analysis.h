/*
	analysis
		for analysis of monte carlo data, and newton-raphson data
*/

#ifndef __ANALYSIS_H_INCLUDED__
#define __ANALYSIS_H_INCLUDED__

#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h> 	// random number generators
#include "simple.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. MonteCarloData class
	2. NewtonRaphsonDatum class
	3. NewtonRaphsonData class
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. MonteCarloData class
-------------------------------------------------------------------------------------------------------------------------*/

class MonteCarloData {
public:
	// constructors, destructor
	MonteCarloData(const string&);
	~MonteCarloData();
	
	// save and load binary
	void loadData(const string&);
	void saveData(const string&) const;
	void saveCorrelator(const string&) const;
	void saveCorrelatorAppend(const string&) const;
	
	// save and load ascii
	void saveResults(const string&) const;
	void loadDataAscii(const string&);
	void saveDataAscii(const string&) const;
	void saveCorrelatorAscii(const string&) const;
	void saveCorrelatorAppendAscii(const string&) const;
	
	// calculations
	void					calcMeans(number& mean, number& meanSqrd);
	void					calcMeans();
	void					calcCorrs(vector<number>& correlator, number& intCorrTime, \
									number& expCorrTime, number& corrErrorSqrd, uint start=0, uint end=0);
	void					calcCorrs(vector<number>& correlator, uint start=1, uint end=0);
	void					calcCorrs(number& intCorrTime, number& expCorrTime, number& corrErrorSqrd, uint start=0, uint end=0);
	void					calcCorrs(uint start=1, uint end=1);
	number					calcJacknife();
	number					calcBootstrap(uint N, const uint& Seed);
	
private:
	uint 			Size;
	gsl_rng* 		Generator;
	number 			Mean;
	number 			MeanSqrd;
	number 			IntCorrTime;
	number 			ExpCorrTime;
	number 			CorrErrorSqrd;
	number 			Jacknife;
	number 			Bootstrap;
	vector<number> 		DataArray;
	vector<number> 		Correlator;
	void			zeroNums();
};

/*-------------------------------------------------------------------------------------------------------------------------
	2. NewtonRaphsonDatum class
-------------------------------------------------------------------------------------------------------------------------*/

class NewtonRaphsonDatum {
public:
	// constructors, destructor
	NewtonRaphsonDatum(const uint& idsize, const uint& datumsize);
	NewtonRaphsonDatum(const vector<string>&, const Parameters&, const vector<number>&);
	~NewtonRaphsonDatum();
	
	// save and load ascii
	void save(const string&) const;
	void load(const string&);
	void load(const vector<string>&);
	void load(const vector<string>&, const Parameters&, const vector<number>&);
	
	// get
	uint size() const;
	vector<string> id() const;
	Parameters parameters() const;
	vector<number> datum() const;
	vector<string> strings() const;
	
	// checks
	bool checkID(const vector<string>&) const;
	bool checkParameters(const Parameters&) const;
	
private:
	uint IDSize;
	vector<string> ID;
	Parameters P;
	uint DatumSize;
	vector<number> Datum;
	
	vector<string> stringVector() const;
};

// operator==
bool operator==(const NewtonRaphsonDatum& lhs, const NewtonRaphsonDatum& rhs);

// operator<<
ostream& operator<<(ostream&, const NewtonRaphsonDatum&);

// operator>>
istream& operator>>(istream&, NewtonRaphsonDatum&);

/*-------------------------------------------------------------------------------------------------------------------------
	3. NewtonRaphsonData class
		
-------------------------------------------------------------------------------------------------------------------------*/

class NewtonRaphsonData {
public:
	// constructors, destructor
	NewtonRaphsonData(const string&, const uint& idsize, const uint& datumsize);
	NewtonRaphsonData(const vector<NewtonRaphsonDatum>&);
	~NewtonRaphsonData();
	
	// save and load ascii
	void save(const string&) const;
	void saveAppend(const string&) const;
	void load(const string&, const uint& idsize, const uint& datumsize);
	
	// get
	uint size() const;
	
	// search
	bool find(const vector<string>& id) const;
	bool find(const Parameters&) const;
	bool find(const vector<string>& id, const Parameters&) const;
	bool find(NewtonRaphsonDatum&) const;
	
private:
	uint Size;
	vector<NewtonRaphsonDatum> DataArray;
};

#endif // __ANALYSIS_H_INCLUDED__
