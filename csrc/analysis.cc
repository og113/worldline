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
	Generator = gsl_rng_alloc(gsl_rng_taus);
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
	saveVectorBinaryAppend<number>(f,Correlator);
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
								 number& expCorrTime, number& corrErrorSqrd) {
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
	for (uint t=1; t<(Size-1); t++) {
		for (uint l=0; l<(Size-t); l++)
			Correlator[t] += DataArray[l]*DataArray[l+t];
		Correlator[t] /= (number)(Size-t);
		Correlator[t] -= Mean*Mean;
		Correlator[t] /= scaling;
		
		IntCorrTime += Correlator[t];
		
		if (Correlator[t]>0 && Correlator[t]<Correlator[t-1] && expBool) {
			expCount++;
			ExpCorrTime += -(number)t/gsl_sf_log(Correlator[t]);
		}
		else
			expBool = false;
	}
	if (expCount>1)
		ExpCorrTime /= (number)expCount;
		
	CorrErrorSqrd = 2.0*IntCorrTime*scaling/(number)(Size-1.0);
	
	// assigning results
	intCorrTime = IntCorrTime;
	expCorrTime = ExpCorrTime;
	correlator = Correlator;
	corrErrorSqrd = CorrErrorSqrd;
	
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs(vector<number>& correlator) {
	number dross1, dross2, dross3;
	calcCorrs(correlator,dross1,dross2,dross3);
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs(number& intCorrTime, number& expCorrTime, number& corrErrorSqrd) {
	vector<number> drossVector;
	calcCorrs(drossVector,intCorrTime,expCorrTime,corrErrorSqrd);
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
	gsl_rng_set(Generator,Seed);
	number mean = 0.0, meanSqrd = 0.0;
	number mean_local, meanSqrd_local;
	uint bootstraps = N*Size;
	for (uint j=0; j<bootstraps; j++) {
		mean_local = 0.0;
		meanSqrd_local = 0.0;
		for (uint j=0; j<Size; j++) {
			loc = (uint)(gsl_rng_uniform (Generator)*N);
			mean_local += DataArray[loc];
			meanSqrd_local += DataArray[loc]*DataArray[loc];
		}
		mean += mean_local;
		meanSqrd += meanSqrd_local;
	}
	mean /= (number)(Size*bootstraps);
	meanSqrd /= (number)(Size*bootstraps);
	return sqrt(meanSqrd-mean*mean);	
}


