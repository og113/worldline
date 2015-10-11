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

// save - saveVectorBinary
template <class T>
void saveVecBinary(const string& f,  const vector<T>& v);

/*-------------------------------------------------------------------------------------------------------------------------
	2. load
-------------------------------------------------------------------------------------------------------------------------*/

// loadVectorBinary
template <class T>
void loadVectorBinary(const string& f, vector<T>& v);


#endif // __PRINT_H_INCLUDED__
