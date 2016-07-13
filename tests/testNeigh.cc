
#include <iostream>
#include <iomanip>
#include <vector>
#include "simple.h"
#include "genloop.h"

using namespace std;

int main() {

cout << "testNeigh: " << endl;

# define dim 4
uint N = 8;

cout << "dim = " << dim << endl;
cout << "N = " << N << endl << endl;

cout << left;
cout << setw(20) << "j" << setw(20) << "posNeigh" << setw(20) << "posNeighDisjoint" << endl;
for (uint j=0; j<N; j++) {
	cout << setw(20) << j << setw(20) << posNeigh(j,N) << setw(20) << posNeighDisjoint(j,N) << endl;
}

cout << endl;

return 0;
}
