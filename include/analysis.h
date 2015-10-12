/*
	analysis
		for analysis of monte carlo data
*/

#ifndef __ANALYSIS_H_INCLUDED__
#define __ANALYSIS_H_INCLUDED__

#include <iostream>
#include <string>
#include <vector>
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
	
	// save and load
	void loadData(const string&);
	void saveData(const string&) const;
	void saveCorrelator(const string&) const;
	void saveResults(const string&) const;
	
	// calculations
	void					calcMeans(number& mean, number& meanSqrd);
	void					calcMeans();
	void					calcCorrs(vector<number>& correlator, number& intCorrTime, \
									number& expCorrTime, number& corrErrorSqrd);
	void					calcCorrs(vector<number>& correlator);
	void					calcCorrs(number& intCorrTime, number& expCorrTime, number& corrErrorSqrd);
	void					calcCorrs();
	number					calcJacknife();
	number					calcBootstrap(const uint& N, const uint& Seed);
	
private:
	uint Size;
	number Mean;
	number MeanSqrd;
	number IntCorrTime;
	number ExpCorrTime;
	number CorrErrorSqrd;
	number Jacknife;
	number Bootstrap;
	vector<number> DataArray;
	vector<number> Correlator;
	void		zeroNums();
};
#endif // __ANALYSIS_H_INCLUDED__
