
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "simple.h"
#include "print.h"
#include "genloop.h"

using namespace std;

int main() {

cout << "testDistances" << endl;

#define dim 4

Point<dim> p1, p2;
p2[0] = 1.5;
p2[dim-1] = 1.5;

cout << "p1 = " << p1 << endl;

cout << "p2 = " << p2 << endl;
	
cout << "Distance = " << Distance(p1,p2) << endl;
cout << "SpatialDistance = " << SpatialDistance(p1,p2) << endl;
cout << "ThermalDistance(beta=1) = " << ThermalDistance(p1,p2,1.0) << endl;

return 0;
}
