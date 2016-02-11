
#include <iostream>
#include <iomanip>
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


cout << setw(12) << "pos" << setw(12) << "npos" << setw(12) << "test" << endl;

int test, pos = 0, Nratio, local, npos, k;
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
	
	k=0;
	npos=pos;
	Nratio = 1;
	while (k<D && npos==pos) {
		if (pos%(Nratio*Nvec[k])!=0)
			npos -= Nratio;
		Nratio *= Nvec[k];
		k++;
	}
	
	cout << setw(12) << pos << setw(12) << npos << setw(12) << test << endl;
	
	pos++;
}

return 0;
}
