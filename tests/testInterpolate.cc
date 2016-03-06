/*
	testInterpolate
		visual test of loop interpolate function
*/


#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include "genloop.h"
#include "folder.h"
#include "simple.h"
#include "nrloop.h"
#include "print.h"

using namespace std;

//dimension
#define dim 4

int main() {

Filename f5 = "data/circle/loops/dim_4/K_5/loop_R_1_rank_0.dat";

Filename fo4 = "data/temp/f4.dat"; 
Filename fo5 = "data/temp/f5.dat";
Filename fo6 = "data/temp/f6.dat";

cout << endl << "testInterpolate: " << endl;
cout << endl << "loading from: " << endl << f5 << endl;
cout << endl << "saving to: " << endl << fo4 << endl << fo5 << endl << fo6 << endl << endl;

Loop<dim> l4(4,0), l5(5,0), l6(6,0);

if (!f5.exists()) {
	cerr << "file not all found: " << endl;
	cerr << f5 << endl ;
	return 1;
}

l5.load(f5);

interpolate(l5,l4);
interpolate(l5,l6);

l4.saveAscii(fo4);
l5.saveAscii(fo5);
l6.saveAscii(fo6);

return 0;
}
