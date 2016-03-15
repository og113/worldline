/*
	testUint
		quick test of unsigned versus signed integer calculations
*/

#include <iostream>
#include "simple.h"

using namespace std;

int main() {

uint u2 = 2, u1 = 1;
number n = -(int)u2/(number)u1;
uint v = 3;
int j = 1;

cout << u1 << endl;
cout << u2 << endl;
cout << j << endl;
cout << n << endl;
cout << v << endl;
cout << (int)(n*j) << endl;
v += (int)(n*j);
cout << v << endl;
v = 3;
n *= -1;
j *= -1;
v += (int)n*j;
cout << (int)n*j << endl;
cout << v << endl;
v = 3;
n *= -1;
v += n*u1;
cout << n*u1 << endl;
cout << v << endl;

return 0;
}
