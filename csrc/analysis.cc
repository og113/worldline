/*
	analysis
		for analysis of monte carlo data
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <gsl/gsl_sf_log.h>
#include "analysis.h"
#include "print.h"
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
	Mean(), MeanSqrd(), IntCorrTime(), ExpCorrTime(), Jacknife(), Bootstrap() {
	loadVectorBinary<number>(f,DataArray);
	Size = DataArray.size();
	Correlator.resize(Size,0.0);
}

// destructor
MonteCarloData::~MonteCarloData() {}

// zero
void MonteCarloData::zeroNums() {
	Mean = 0.0;
	MeanSqrd = 0.0;
	IntCorrTime = 0.0;
	ExpCorrTime = 0.0;
	Jacknife = 0.0;
	Bootstrap = 0.0;
}

// load Data
void MonteCarloData::loadData(const string& f) {
	loadVectorBinary<number>(f,DataArray);
	Size = DataArray.size();
	Correlator.resize(Size,0.0);
	zeroNums();
}

// save Data
void MonteCarloData::saveData(const string& f) const {
	saveVectorBinary<number>(f,DataArray);
}

// save Correlator
void MonteCarloData::saveCorrelator(const string& f) const {
	if (abs(IntCorrTime)<MIN_NUMBER && abs(ExpCorrTime)<MIN_NUMBER) {
		cerr << "MonteCarloData::saveCorrelator error: correlator not calculated yet" << endl;
		return;
	}
	saveVectorBinary<number>(f,Correlator);
}

// save Results
void MonteCarloData::saveResults(const string& f) const {
	if (abs(Mean)<MIN_NUMBER && abs(MeanSqrd)<MIN_NUMBER) {
		cerr << "MonteCarloData::saveResults error: results not calculated yet" << endl;
		return;
	}
	FILE * os;
	os = fopen(f.c_str(),"a");
	fprintf(os,"%12i%24g%24g",Size,Mean,MeanSqrd);
	if (abs(IntCorrTime)>MIN_NUMBER || abs(ExpCorrTime)>MIN_NUMBER)
		fprintf(os,"%24g%24g",IntCorrTime,ExpCorrTime);
	if (abs(Jacknife)>MIN_NUMBER)
		fprintf(os,"%24g",Jacknife);
	if (abs(Bootstrap)>MIN_NUMBER)
		fprintf(os,"%24g",Bootstrap);
	fprintf(os,"\n");
	fclose(os);
}

// calc Mean and MeanSqrd
void MonteCarloData::calcMeans(number& mean, number& meanSqrd) {
	if (Size==0) {
		cerr << "MonteCarloData::calcMeans error: no data" << endl;
	}
	if (abs(Mean)>MIN_NUMBER || abs(MeanSqrd)>MIN_NUMBER) {
		mean = Mean;
		meanSqrd = MeanSqrd;
		return;
	}
	for (uint j=0; j<Size; j++) {
		Mean += DataArray[j];
		MeanSqrd += DataArray[j]*DataArray[j];
	}
	Mean /= (number)Size;
	MeanSqrd /= (number)Size;
	mean = Mean;
	meanSqrd = MeanSqrd;
}

// calc Mean and MeanSqrd
void MonteCarloData::calcMeans() {
	number dross1, dross2;
	calcMeans(dross1,dross2);
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs(vector<number>& correlator, number& intCorrTime,\
								 number& expCorrTime, number& corrError) {
	if (Size==0) {
		cerr << "MonteCarloData::calcCorrs error: no data" << endl;
	}
	if (abs(IntCorrTime)>MIN_NUMBER || abs(ExpCorrTime)>MIN_NUMBER) {
		intCorrTime = IntCorrTime;
		expCorrTime = ExpCorrTime;
		correlator = Correlator;
		return;
	}
	if (abs(Mean)<MIN_NUMBER && abs(MeanSqrd)<MIN_NUMBER) {
		calcMeans();
		return;
	}
	// scaling for correlator
	number scaling = MeanSqrd-Mean*Mean;
	Correlator[0] = 1.0;
	
	uint expCount = 0;
	bool expBool = true;
	
	IntCorrTime = 0.5;
	// double sum to calculate correlator
	for (uint k=1; k<(Size-1); k++) {
		for (uint l=0; l<(Size-k); l++)
			Correlator[k] += DataArray[l]*DataArray[l+k];
		Correlator[k] /= (number)(Size-1.0-k);
		Correlator[k] -= Mean*Mean;
		Correlator[k] /= scaling;
		
		IntCorrTime += Correlator[k];
		
		if (Correlator[k]>0 && Correlator[k]<Correlator[k-1] && expBool) {
			expCount++;
			ExpCorrTime += -(number)k/gsl_sf_log(Correlator[k]);
		}
		else
			expBool = false;
	}
	if (expCount!=0)
		ExpCorrTime /= (number)expCount;
		
	CorrError = 2.0*IntCorrTime*scaling/(number)(Size-1.0);
	CorrError = sqrt(CorrError);
	
	// assigning results
	intCorrTime = IntCorrTime;
	expCorrTime = ExpCorrTime;
	correlator = Correlator;
	corrError = CorrError;
	
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs(vector<number>& correlator) {
	number dross1, dross2, dross3;
	calcCorrs(correlator,dross1,dross2,dross3);
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs() {
	number dross1, dross2, dross3;
	vector<number> drossVector;
	calcCorrs(drossVector,dross1,dross2,dross3);
}

// calcJacknife
number MonteCarloData::calcJacknife() {
	if (abs(Mean)<MIN_NUMBER && abs(MeanSqrd)<MIN_NUMBER) {
		calcMeans();
	}
	if (abs(IntCorrTime)<MIN_NUMBER && abs(ExpCorrTime)<MIN_NUMBER) {
		calcCorrs();
	}
	cerr << "MonteCarloData::Jacknife error: haven't written function yet" << endl;
	return 0.0;	
}

// calcBootStrap
number MonteCarloData::calcBootstrap(const uint& N, const uint& Seed) {
	if (abs(Mean)<MIN_NUMBER && abs(MeanSqrd)<MIN_NUMBER) {
		calcMeans();
	}
	cerr << "MonteCarloData::Bootstrap error: haven't written function yet" << endl;
	return 0.0;	
}


