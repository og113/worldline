/*
	analysis
		for analysis of monte carlo data
*/

#include <iostream>
#include <string>
#include <vector>
#include "analysis.h"
#include "simple.h"

/*-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. MonteCarloData class
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. MonteCarloData class
		
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
MonteCarloData::MonteCarloData(const string& f): 
	Mean(), MeanSqrd(), IntCorrTime(), ExpCorrTime(), Jacknife(), BootStrap() {
	loadVectorBinary<number>(f,DataArray);
	Size = DataArray.size();
	Correlator.resize(Size,0.0);
}

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
	number Jacknife;
	number BootStrap;
	vector<number> DataArray;
	vector<number> Correlator;
};

