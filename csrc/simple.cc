/*-------------------------------------------------------------------------------------------------------------------------
definitions of some very simple functions and classes
-------------------------------------------------------------------------------------------------------------------------*/

#include <cstring>
#include <cstdlib> //for rand, srand
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <sys/time.h>
#include <glob.h>
#include "simple.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - number to string and string to number
	2 - absDiff
	3 - factorial
	4 - currentDateTime, currentPartSec
	5 - copyFile
	6 - count in files
	7 - smallestLoc
	8 - randDouble
	9 - mod
	10 - delta
	11 - explicit instantiation
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. number to string and string to number
		- number to string
		- string to number
		- isNumber
-------------------------------------------------------------------------------------------------------------------------*/

//to convert number to string, usage is string str = NumberToString<number type>(x);
template <class T>
string numberToString ( const T& Number ) {
	stringstream ss;
	ss << Number;
	return ss.str();
}
	
//shorthand version;
template <class T>
string nts ( const T& Number ) {
	stringstream ss;
	ss << Number;
	return ss.str();
}

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <class T>
T stringToNumber ( const string& Text ) {
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}
	
//shorthand version;
template <class T>
T stn ( const string& Text ) {
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

// is number?
bool isNumber( const string& Text ) {
	return( strspn( Text.c_str(), "0123456789" ) == Text.size() );
}
	
/*-------------------------------------------------------------------------------------------------------------------------
	2. absDiff
		- double
		- comp
		- vec
		- cVec
		- mat
		- cMat
-------------------------------------------------------------------------------------------------------------------------*/
	
// absDiff double
double absDiff(const double& numA, const double& numB) {
	if (abs(numA)>MIN_NUMBER && abs(numB)>MIN_NUMBER) return 2.0*abs(numA-numB)/sqrt(numA*numA+numB*numB);
	else if (abs(numA)>MIN_NUMBER || abs(numB)>MIN_NUMBER) return abs(numA-numB);
	else return 0.0;
}

// absDiff comp
double absDiff (const comp& numA, const comp& numB) {
	if (abs(numA)>MIN_NUMBER && abs(numB)>MIN_NUMBER) return 2.0*abs(numA-numB)/sqrt(norm(numA)+norm(numB));
	else if (abs(numA)>MIN_NUMBER || abs(numB)>MIN_NUMBER) return abs(numA-numB);
	else return 0.0;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. factorial
-------------------------------------------------------------------------------------------------------------------------*/

//factorial
int factorial(const int& f_input){
	int f_result = 1;
	for (int l=0; l<f_input; l++) f_result *= (l+1);
	return f_result;
}

/*-------------------------------------------------------------------------------------------------------------------------
	4. currentDateTime, currentPartSec
	
N.B. Visit http://en.cppreference.com/w/cpp/chrono/c/strftime for more information about date/time format
-------------------------------------------------------------------------------------------------------------------------*/

//getting the date and time
string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%y%m%d%H%M%S", &tstruct);
    return buf;
}

//getting the part of the second, a number between 0 and 1e6, divide by 1e6 to get the fraction
string currentPartSec() {
    timeval tim;
    gettimeofday(&tim, NULL);
    return numberToString<double>(tim.tv_usec);
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. copyFile, fileExists
-------------------------------------------------------------------------------------------------------------------------*/

//copy a file
void copyFile(const string & inputFile, const string & outputFile) {
	ifstream  is(inputFile.c_str(), ios::binary);
	ofstream  os(outputFile.c_str(), ios::binary);
	os << is.rdbuf();
}

// fileExists
bool fileExists(const string& f) {
	glob_t globbuf;
	int err = glob(f.c_str(), GLOB_NOSORT, NULL, &globbuf); // GLOB_NOSORT if doing sort elsewhere, otherwise set to zero
	if(err == 0) {	
		if(globbuf.gl_pathc>0)
			globfree(&globbuf);
		return true;
	}
	else if (err==GLOB_NOSPACE) {
        	cerr << "Folder search error " << err << ", running out of memory" << endl;
		return false;
	}
	else if (err==GLOB_ABORTED) {
        	cerr << "Folder search error " << err << ", GLOB_ABORTED" << endl;
		return false;
	}
	else if (err==GLOB_NOMATCH) {
		return false;
	}
	else {
		cerr << "Folder search error " << err << ", unknown" << endl;
		return false;
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	6. countLines, countColumns, countDoubles
-------------------------------------------------------------------------------------------------------------------------*/

// count lines
//count non-empty lines of a file
uint countLines(const string & file_to_count) {
	ifstream fin;
	fin.open(file_to_count.c_str(), ios::in);
	if (!fin.good()) cerr << "countLines error: " << file_to_count << " not opened properly." << endl;
	string line;
	unsigned int counter = 0;
	while(!fin.eof()) {
		getline(fin,line);
		if(line.empty()) continue;
		counter++;
	}	
	fin.close();
    return counter;
}

// countColumns
uint countColumns(const string & file_to_count) {
	ifstream fin;
	fin.open(file_to_count.c_str(), ios::in);
	if (!fin.good()) cerr << "countRows error: " << file_to_count << " not opened properly." << endl;
	string line;
	unsigned int counter = 0;
	while(!fin.eof()) {
		getline(fin,line);
		if(line.empty()) continue;
		else {
			stringstream ss(line);
			string temp;
			while (ss >> temp) {
				counter++;
			}
			break;
		}
	}		
	fin.close();
    return counter;
}

// countDoubles
uint countDoubles(const string& f) {
	uint lines = -1; // for some reason we should start on -1 not 0, see testBinaryPrint for verification
	ifstream is;
	is.open(f.c_str(),ios::binary);
	double dross;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&dross),sizeof(double));
		lines++;
	}
	is.close();
	return lines;
}

// count comp (binary)
uint countComp(const string& f) {
	uint lines = -1; 
	ifstream is;
	is.open(f.c_str(),ios::binary);
	comp dross;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&dross),sizeof(comp));
		lines++;
	}
	is.close();
	return lines;
}

// count type (binary)
template <class T>
uint countType(const string& f, const T& t) {
	uint lines = -1; 
	ifstream is;
	is.open(f.c_str(),ios::binary);
	T dross;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&dross),sizeof(T));
		lines++;
	}
	is.close();
	return lines;
}

/*-------------------------------------------------------------------------------------------------------------------------
	7. smallestLoc
-------------------------------------------------------------------------------------------------------------------------*/

// smallestLoc
// function giving location of smallest element of a vector of type T
template <typename T>
uint smallestLoc(const vector<T>& inVector) {
	uint loc = 0;
	for(uint l=1;l<inVector.size();l++) {
		if (abs(inVector[l])<abs(inVector[loc])) {
			loc = l;
		}
	}
	return loc;
}

/*-------------------------------------------------------------------------------------------------------------------------
	8. randDouble

	n.b. must first seed rand with srand(time(NULL)) or something similar
-------------------------------------------------------------------------------------------------------------------------*/

// randDouble
double randDouble(const double& min, const double& max) {
	double f = (double)rand() / RAND_MAX;
    return min + f*(max - min);
}

/*-------------------------------------------------------------------------------------------------------------------------
	9. mod
-------------------------------------------------------------------------------------------------------------------------*/

// mod
template <class T>
T mod(const T& x, const T& min, const T& max) {
	T Min, Max;
	if (min<max) {
		Min = min;
		Max = max;
	}
	else if (min>max) {
		Min = max;
		Max = min;
	}
	else {
		cerr << "mod error, range of size zero" << endl;
		return 1.0;
	}
		
	if (x>=Min && x<=Max)
		return x;
	else if (x>Max) {
		int ranges = (int)((x-Min)/(Max-Min));
		return x-(T)ranges*(Max-Min);
	}
	else {
		int ranges = (int)((Max-x)/(Max-Min));
		return x+(T)ranges*(Max-Min);
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	10. delta
-------------------------------------------------------------------------------------------------------------------------*/

double delta(const uint& i, const uint& j) {
	return (i==j? 1.0: 0.0);
}

/*-------------------------------------------------------------------------------------------------------------------------
	11. explicit instantiation
		- numberToString, stringToNumber
		- smallestLoc
		- countType
-------------------------------------------------------------------------------------------------------------------------*/

template string numberToString<int>(const int&);
template string numberToString<uint>(const uint&);
template string numberToString<lint>(const lint&);
template string numberToString<long long unsigned>(const long long unsigned&);
template string numberToString<double>(const double&);
template string numberToString<comp>(const comp&);

template string nts<int>(const int&);
template string nts<uint>(const uint&);
template string nts<lint>(const lint&);
template string nts<long long unsigned>(const long long unsigned&);
template string nts<double>(const double&);
template string nts<comp>(const comp&);

template int stringToNumber<int>(const string&);
template uint stringToNumber<uint>(const string&);
template lint stringToNumber<lint>(const string&);
template long long unsigned stringToNumber<long long unsigned>(const string&);
template double stringToNumber<double>(const string&);
template comp stringToNumber<comp>(const string&);

template int stn<int>(const string&);
template uint stn<uint>(const string&);
template lint stn<lint>(const string&);
template long long unsigned stn<long long unsigned>(const string&);
template double stn<double>(const string&);
template comp stn<comp>(const string&);

template uint smallestLoc<int>(const vector<int>&);
template uint smallestLoc<double>(const vector<double>&);
template uint smallestLoc<comp>(const vector<comp>&);

template uint countType<int>(const string&,const int&);
template uint countType<uint>(const string&,const uint&);
template uint countType<double>(const string&,const double&);
template uint countType<comp>(const string&,const comp&);

template double mod<double>(const double& x, const double& min, const double& max);
template int mod<int>(const int& x, const int& min, const int& max);
