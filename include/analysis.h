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
	2. NewtonRaphsonDatum and NewtonRaphsonData classes
	
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
	2. NewtonRaphsonDatum and NewtonRaphsonDatum classes
-------------------------------------------------------------------------------------------------------------------------*/
/*
// operator==
bool operator==(const NewtonRaphsonDatum& lhs, const NewtonRaphsonDatum& rhs);

// operator<<
ostream& operator<<(ostream&, const NewtonRaphsonDatum&);

class NewtonRaphsonDatum {
public:
	// constructors, destructor
	NewtonRaphsonDatum(const vector<string>&, const Parameters&, const vector<number>&);
	~NewtonRaphsonDatum();
	
	// save and load ascii
	void save(const string&) const;
	void load(const string&);
	
	// get
	vector<string> id() const;
	Parameters parameters() const;
	vecor<number> datum() const;
	
	// checks
	bool checkID(const vector<string>&) const;
	bool checkParameters(const Parameters&) const;
	
	// ==
	//friend bool operator==(const NewtonRaphsonDatum& lhs, const NewtonRaphsonDatum& rhs);
	
	// operator<<
	//friend ostream& operator<<(ostream&, const NewtonRaphsonDatum&);
	
private:
	uint IDSize;
	uint DatumSize;
	vector<string> ID;
	Parameters P;
	vector<number> Datum;
};

class NewtonRaphsonData {
public:
	// constructors, destructor
	NewtonRaphsonData();
	NewtonRaphsonData(const string&);
	NewtonRaphsonData(const vector<NewtonRaphsonDatum>&);
	~NewtonRaphsonData();
	
	// save and load ascii
	void save(const string&) const;
	void load(const string&);
	
	// search
	bool find(const Paramters&);
	bool find(const Paramters&, NewtonRaphsonDatum&);
	
private:
	uint Size;
	vector<NewtonRaphsonDatum> DataArray;
};
*/
#endif // __ANALYSIS_H_INCLUDED__
