/*
	quick program to output suitable ranges for highTemp
*/

#include <iostream>
#include <fstream>
#include <string>
#include "simple.h"
#include "folder.h"
#include "parameters.h"

using namespace std;

int main() {
// loading parameters
ParametersRange pr;
Parameters p;
pr.load("inputs4");
p = pr.Min;

// parameter loops
uint Npl = 1; // number of parameter loops
Parameters::Label label = static_cast<Parameters::Label>(0);
if (pr.toStep(label)) {
	Npl = (pr.Steps)[label-1];
	cout << "looping " << label << " over " << Npl << " steps" << endl;
}

for (uint pl=0; pl<Npl; pl++) {
	// stepping parameters
	if (pr.toStep(label) && pl>0)
		p.step(pr);
	cout << "step " << pl << endl;
	cout << p << endl;
}

}
