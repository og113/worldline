
#include <iostream>
#include <vector>

using namespace std;

int main() {

int D, N = 1;
cout << "D: ";
cin >> D;

vector<int> Nvec(D), nvec(D);

for (int j=0; j<D; j++) {
	cout << "N" << j << ": ";
	cin >> Nvec[j];
	N *= Nvec[j];
}


int test, pos = 0, Nratio, local;
for (int j=0; j<N; j++) {
	test = -j;
	local = pos;
	Nratio = N;
	for (int l=0; l<D; l++) {
		Nratio /= Nvec[D-l-1];
		nvec[D-l-1] = local/Nratio;
		local -= nvec[D-l-1]*Nratio;
		test += nvec[D-l-1]*Nratio;
	}
	pos++;
	cout << "test " << j << ": " << test << endl;
}

return 0;
}
