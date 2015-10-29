/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
  
#ifndef __PRINT_H_INCLUDED__
#define __PRINT_H_INCLUDED__

#include "simple.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. save
	2. load
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. save
-------------------------------------------------------------------------------------------------------------------------*/

// save - saveVectorAscii
template <class T>
void saveVectorAscii(const string& f,  const vector<T>& v);

// save - saveVectorAsciiAppend
template <class T>
void saveVectorAsciiAppend(const string& f,  const vector<T>& v);

// save - saveVectorBinary
template <class T>
void saveVectorBinary(const string& f,  const vector<T>& v);

// save - saveVectorBinaryAppend
template <class T>
void saveVectorBinaryAppend(const string& f,  const vector<T>& v);

/*-------------------------------------------------------------------------------------------------------------------------
	2. load
-------------------------------------------------------------------------------------------------------------------------*/

// loadVectorBinary
template <class T>
void loadVectorBinary(const string& f, vector<T>& v);

// loadVectorAscii
template <class T>
void loadVectorAscii(const string& f, vector<T>& v);

// loadVectorAsciiColumn
template <class T>
void loadVectorAsciiColumn(const string& f, vector<T>& v, const uint& col);

#endif // __PRINT_H_INCLUDED__
