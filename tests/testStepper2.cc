/*
	main for program to test stepper.cc
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "stepper.h"
#include "simple.h"
#include "print.h"
#include "folder.h"

using namespace std;

double f(const Point2d& p) {
	//return pow(p.X,2.0)+pow(p.Y,2.0);
	//return pow(p.X,-2.0)+pow(p.Y,-2.0);
	//return p.X*exp(-pow(p.X,2.0)-pow(p.Y,2.0));
	return exp(+pow(p.X,2.0)+pow(p.Y,2.0))/p.X/p.Y;
}


int main() {
cout << "test stepper2: " << endl;

uint loops = 20;
bool coutEveryLoop = true;
bool coutAngle = true;

string stoFile = "step0";
StepperOptions sto;
ifstream is;
is.open(stoFile.c_str());
if (is.good()) {
	is >> sto;
}
else {
	cerr << "error: cannot open " << stoFile << endl;
}

Point2d P(1.0,1.0);
double maxError = 0.0;
double f0 = f(P);

Stepper st(sto,P);
st.addResult(f(P));

string dataFile = "tests/data/testStepper.dat";
ofstream os;
os.open(dataFile.c_str());

cout << "StepperOptions:" << endl << sto << endl;

if (coutEveryLoop) {
	cout << P << setw(15) << 0.0 << setw(15) << f(P) << setw(15) << "y" << endl;
}
for (unsigned int j=0; j<loops; j++) {
	if (st.keep()) {
		os << P << setw(15) << f(P) << endl;
		if (!coutEveryLoop && coutAngle)
			cout << st.stepAngle() << endl;
		maxError = (absDiff(f(P),f0)>maxError? absDiff(f(P),f0): maxError);
	}
	st.step();
	if (coutEveryLoop) {
		string keep = (st.keep()? "y": "n");
		cout << P << setw(15) << st.stepAngle() << setw(15) << f(P);
	}
	P = st.point();
	st.addResult(f(P));
	if (coutEveryLoop) {
	string keep = (st.keep()? "y": "n");
	cout << setw(15) << keep << endl;
	}
}
os.close();

string stepperFile = "tests/data/stepper.txt";
string stepperFile2 = "tests/data/stepper2.txt";
st.save(stepperFile);
cout << "saved stepper in " << stepperFile << endl;

Stepper st2(sto,P);
st2.load(stepperFile);
cout << "loaded stepper from " << stepperFile << endl;
st2.save(stepperFile2);
cout << "saved stepper in " << stepperFile2 << endl;

return 0;
}
