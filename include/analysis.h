/*
	analysis
		for analysis of monte carlo data
*/

#ifndef __ANALYSIS_H_INCLUDED__
#define __ANALYSIS_H_INCLUDED__

#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h> 	// random number generators
#include "simple.h"

/*-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. MonteCarloData class
	
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
#endif // __ANALYSIS_H_INCLUDED__
