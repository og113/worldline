/*-------------------------------------------------------------------------------------------------------------------------
declarations of some very simple functions and classes
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __SIMPLE_H_INCLUDED__
#define __SIMPLE_H_INCLUDED__

#include <string>
#include <vector>
#include <complex>

using namespace std;

#ifndef PI
#define PI 3.1415926535897932
#endif
#ifndef SQRT2
#define SQRT2 1.4142135623730950
#endif
#ifndef MIN_NUMBER
#define MIN_NUMBER 1.0e-16
#endif
#ifndef LOG_MIN_NUMBER
#define LOG_MIN_NUMBER -36.841361488 // calculated as ln(1.0e-16)
#endif

typedef double number;
typedef unsigned int uint;
typedef unsigned long int lint;
typedef complex<double> comp;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - number to string and string to number
	2 - absDiff
	3 - factorial
	4 - currentDateTime, currentPartSec
	5 - copyFile, fileExists
	6 - count in files
	7 - smallestLoc
	8 - randDouble
	9 - mod
	10 - delta
	11 - sigFig
	12 - sign
	13 - split
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. number to string and string to number
		- number to string
		- string to number
		- isNumber
		- string to bool
-------------------------------------------------------------------------------------------------------------------------*/

template <class T>
string numberToString ( const T& Number );

template <class T>
string nts ( const T& Number );

template <class T>
string nts ( const T& Number, const uint& prec);

template <class T>
T stringToNumber ( const string& Text );

template <class T>
T stn ( const string& Text );

bool isNumber( const string& Text );

/*-------------------------------------------------------------------------------------------------------------------------
	2. absDiff
		- double
		- comp
		- vec
		- cVec
		- mat
		- cMat
-------------------------------------------------------------------------------------------------------------------------*/

double absDiff (const double& A, const double& B);
double absDiff (const comp& A, const comp& B);

/*-------------------------------------------------------------------------------------------------------------------------
	3. factorial
-------------------------------------------------------------------------------------------------------------------------*/

//factorial
int factorial(const int& f_input);

/*-------------------------------------------------------------------------------------------------------------------------
	4. currentDateTime, currentPartSec
-------------------------------------------------------------------------------------------------------------------------*/

//getting the date and time
string currentDateTime();

// getting the part of the second
string currentPartSec();

/*-------------------------------------------------------------------------------------------------------------------------
	5. copyFile
-------------------------------------------------------------------------------------------------------------------------*/

//copy a file
void copyFile(const string & inputFile, const string & outputFile);

// fileExists
bool fileExists(const string& f);

/*-------------------------------------------------------------------------------------------------------------------------
	6. count in file
		- countLines
		- countColumns
		- countDoubles
		- countComp
		- countType
-------------------------------------------------------------------------------------------------------------------------*/

// count lines
uint countLines(const string & f);

// count columns
uint countColumns(const string& f);

// count columns
uint countColumns(const string& f, const char& separator);

// count doubles (binary)
uint countDoubles(const string& f);

// count comp (binary)
uint countComp(const string& f);

// count type (binary)
template <class T>
uint countType(const string& f, const T& t);

/*-------------------------------------------------------------------------------------------------------------------------
	7. smallestLoc
-------------------------------------------------------------------------------------------------------------------------*/

// smallestLoc
template <class T>
uint smallestLoc(const vector<T>&);

/*-------------------------------------------------------------------------------------------------------------------------
	8. randDouble
	
	n.b. must first seed rand with srand(time(NULL)) or something similar
-------------------------------------------------------------------------------------------------------------------------*/

// randDouble
double randDouble(const double& min, const double& max);

/*-------------------------------------------------------------------------------------------------------------------------
	9. mod
-------------------------------------------------------------------------------------------------------------------------*/

// mod template
template <class T>
T mod(const T& x, const T& min, const T& max);

/*-------------------------------------------------------------------------------------------------------------------------
	10. delta
-------------------------------------------------------------------------------------------------------------------------*/

double delta(const uint&, const uint&);

/*-------------------------------------------------------------------------------------------------------------------------
	11. sigFig
-------------------------------------------------------------------------------------------------------------------------*/

double sigFig(double num, double N);

/*-------------------------------------------------------------------------------------------------------------------------
	12. sign
-------------------------------------------------------------------------------------------------------------------------*/

double sign(double num);

/*-------------------------------------------------------------------------------------------------------------------------
	13. splitString
-------------------------------------------------------------------------------------------------------------------------*/

vector<string> splitString(const string & str, const string & delimiters);

#endif // __SIMPLE_H_INCLUDED__
