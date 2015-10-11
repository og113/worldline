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
	void save(const string&) const;
	void load(const string&);
	
	// calculations
	void					calcMeans(number& mean, number& meanSqrd);
	void					calcCorrs(vector<number>& correlator, const number& intCorrTime, const number& expCorrTime);
	number					calcJacknife();
	number					calcBootsrap(const uint&);
	
private:
	number Size;
	number Mean;
	number MeanSqrd;
	number IntCorrTime;
	number ExpCorrTime;
	number ErrorJacknife;
	number ErrorBootStrap;
	vector<number> DataArray;
	vector<number> Correlator;
};
#endif // __ANALYSIS_H_INCLUDED__
