/*
	quick test of genFloop
*/


#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>       /* time */
#include <gsl/gsl_sf_exp.h>
#include "genfloop.h"
#include "simple.h"


using namespace std;

int main() {

#define dim 2

uint K = 2, Seed = time(NULL)+3;

FLoop<2> fl(K,Seed);

fl.grow();

cout << S0(fl) << endl;


return 0;
}
