/*
	perturbative finite temperature action
*/

#include <cmath>
#include <iostream>
#include <string>
#include "simple.h"

using namespace std;

int main(int argc, char** argv) {

number x, T, result = PI-1.0/4.0, jx;
uint N = 100;

if (argc<2)  {
	cout << "Type argument, T: ";
	cin >> T;
}
else {
	T = stn<number>(argv[1]);
	if (argc==3)
		N = stn<uint>(argv[2]);
}

if (T>MIN_NUMBER) {
	x = 1.0/T;

	for (uint j=1; j<=N; j++) {
		jx = x*(number)j;
		result -= 0.5*( (-2.0/pow(jx,2) + 1.0)/sqrt(1.0 - 4.0/pow(jx,2)) - 1.0);
	}
}

cout << result << endl;

return 0;
}
