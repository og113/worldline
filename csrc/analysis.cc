/*
	analysis
		for analysis of monte carlo data
*/

#include <iostream>
#include <iomanip>
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
	2. NewtonRaphsonDatum and NewtonRaphsonDatum classes
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. MonteCarloData class
		
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
MonteCarloData::MonteCarloData(const string& f): 
	Mean(), MeanSqrd(), IntCorrTime(), ExpCorrTime(), Jacknife(), Bootstrap() {
	loadVectorBinary(f,DataArray);
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
	loadVectorBinary(f,DataArray);
	Size = DataArray.size();
	Correlator.resize(Size,0.0);
	zeroNums();
}

// load data ascii
void MonteCarloData::loadDataAscii(const string& f) {
	loadVectorAscii(f,DataArray);
}

// save Data
void MonteCarloData::saveData(const string& f) const {
	saveVectorBinary(f,DataArray);
}

// save Correlator
void MonteCarloData::saveCorrelator(const string& f) const {
	if (abs(IntCorrTime)<MIN_NUMBER && abs(ExpCorrTime)<MIN_NUMBER) {
		cerr << "MonteCarloData::saveCorrelator error: correlator not calculated yet" << endl;
		return;
	}
	saveVectorBinary(f,Correlator);
}

// save Correlator Append
void MonteCarloData::saveCorrelatorAppend(const string& f) const {
	if (abs(IntCorrTime)<MIN_NUMBER && abs(ExpCorrTime)<MIN_NUMBER) {
		cerr << "MonteCarloData::saveCorrelator error: correlator not calculated yet" << endl;
		return;
	}
	saveVectorBinaryAppend(f,Correlator);
}

// save Data ascii
void MonteCarloData::saveDataAscii(const string& f) const {
	saveVectorAscii(f,DataArray);
}

// save Correlator ascii
void MonteCarloData::saveCorrelatorAscii(const string& f) const {
	if (abs(IntCorrTime)<MIN_NUMBER && abs(ExpCorrTime)<MIN_NUMBER) {
		cerr << "MonteCarloData::saveCorrelator error: correlator not calculated yet" << endl;
		return;
	}
	saveVectorAscii(f,Correlator);
}

// save Correlator append ascii
void MonteCarloData::saveCorrelatorAppendAscii(const string& f) const {
	if (abs(IntCorrTime)<MIN_NUMBER && abs(ExpCorrTime)<MIN_NUMBER) {
		cerr << "MonteCarloData::saveCorrelator error: correlator not calculated yet" << endl;
		return;
	}
	saveVectorAsciiAppend(f,Correlator);
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
								 number& expCorrTime, number& corrErrorSqrd, uint start, uint end) {
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
	
	if (end==0)
		end = Size-1;
	
	IntCorrTime = 0.5;
	// double sum to calculate correlator
	for (uint t=1; t<(Size-1); t++) {
		for (uint l=0; l<(Size-t); l++)
			Correlator[t] += DataArray[l]*DataArray[l+t];
		Correlator[t] /= (number)(Size-t);
		Correlator[t] -= Mean*Mean;
		Correlator[t] /= scaling;
		
		if (t>start && t<end)
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
		
	CorrErrorSqrd = 2.0*IntCorrTime*scaling/(number)(end-start);
	
	// assigning results
	intCorrTime = IntCorrTime;
	expCorrTime = ExpCorrTime;
	correlator = Correlator;
	corrErrorSqrd = CorrErrorSqrd;
	
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs(vector<number>& correlator, uint start, uint end) {
	number dross1, dross2, dross3;
	calcCorrs(correlator,dross1,dross2,dross3);
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs(number& intCorrTime, number& expCorrTime, number& corrErrorSqrd, uint start, uint end) {
	vector<number> drossVector;
	calcCorrs(drossVector,intCorrTime,expCorrTime,corrErrorSqrd);
}

// calc Correlator and correlation times
void MonteCarloData::calcCorrs(uint start, uint end) {
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

// calcBootStrap - n.b. return square of error
number MonteCarloData::calcBootstrap(uint N, const uint& Seed) {
	if (abs(Mean)<MIN_NUMBER && abs(MeanSqrd)<MIN_NUMBER) {
		calcMeans();
	}
	if (N==0)
		N++;
	gsl_rng_set(Generator,Seed);
	number mean = 0.0, meanSqrd = 0.0;
	number mean_local, meanSqrd_local;
	uint loc;
	for (uint j=0; j<N; j++) {
		mean_local = 0.0;
		meanSqrd_local = 0.0;
		for (uint j=0; j<Size; j++) {
			loc = (uint)(gsl_rng_uniform (Generator)*Size);
			mean_local += DataArray[loc];
			meanSqrd_local += DataArray[loc]*DataArray[loc];
		}
		mean += mean_local;
		meanSqrd += meanSqrd_local;
	}
	mean /= (number)(Size*N);
	meanSqrd /= (number)(Size*N);
	return meanSqrd-mean*mean;	
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. NewtonRaphsonDatum and NewtonRaphsonDatum classes
		
-------------------------------------------------------------------------------------------------------------------------*/
/*
// checks
bool checkID(const vector<string>& id) const {
	for (uint j=0; j<IDSize; j++) {
		if (!(id[j].empty())) {
			if ((id[j]).compare(ID[j])!=0)
				return false;
		}
	}
	return true;
}

// checkParameters
bool checkParameters(const Parameters& p) const {
	return P==p;
}

// ==
bool operator==(const NewtonRaphsonDatum& lhs, const NewtonRaphsonDatum& rhs) {
	return (lhs.checkID(rhs.id()) && lhs.checkParameters(rhs.parameters()));
}

// constructor
NewtonRaphsonDatum::NewtonRaphsonDatum(const string& id, const Parameters& p, const vector<number>& d): 
	ID(id), P(p), Datum(d) {
	IDSize = ID.size();
	DatumSize = Datum.size();
}

// destructor
NewtonRaphsonDatum::~NewtonRaphsonDatum() {}

// save
void NewtonRaphsonDatum::save(const string& f) const {
	uint totalSize = IDSize + Parameters::Size + DatumSize;
	vector<string> PVector = P.valueVector();
	vector<string> toSave(totalSize);
	
	for (uint j=0; j<IDSize; j++)
		toSave[j] = ID[j]; 
	for (uint j=0; j<Parameters::Size; j++)
		toSave[IDSize+j] = PVector[j]; 
	for (uint j=0; j<DatumSize; j++)
		toSave[IDSize+Parameters::Size+j] = nts(Datum[j],16);
		
	saveVectorCsvAppend(f,toSave);
}

// load
void NewtonRaphsonDatum::load(const string& f) {
	uint totalSize = IDSize + Parameters::Size + DatumSize;
	vector<string> totalVector, pVector;
	loadVectorCsvAppend(f,totalVector);
	pVector = vector(totalVector.begin()+1,totalVector.begin()+1+Parameters::Size);
	ID = totalVector[0];
	P.load(pVector);
	for (uint j=0; j<DatumSize; j++)
		Datum[j] = stn(totalVector[IDSize+Parameters::Size+j]);
}

// parameters
Parameters NewtonRaphsonDatum::parameters() const {
	return P;
}

// operator<<
ostream& operator<<(ostream& os, const NewtonRaphsonDatum& d);


// constructor
NewtonRaphsonData::NewtonRaphsonData(): Size(), DatArray() {
}

// constructor
NewtonRaphsonData::NewtonRaphsonData(const string& f): Size(), DatArray() {
	load(f);
}

// constructor
NewtonRaphsonData::NewtonRaphsonData(const vector<NewtonRaphsonDatum>& v) {
	DataArray = v;
	Size = v.size();
}

// destructor
NewtonRaphsonData::~NewtonRaphsonData() {}

// save
void NewtonRaphsonData::save(const string& f) const {
	
}

// load
void NewtonRaphsonData::load(const string& f);

// find
bool find(const Paramters&);

// find
bool find(const Paramters&, NewtonRaphsonDatum&);*/
