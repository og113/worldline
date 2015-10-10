// test of parameters

#include <iostream>
#include <fstream>
#include <string>
#include "parameters.h"

using namespace std;

int main() {

ParametersRange pr;
pr.load("inputs");
pr.save("inputs");
cout << pr << endl;

Parameters p;
p.load("inputs");
cout << p << endl;

return 0;
}
