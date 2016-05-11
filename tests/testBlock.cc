
#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include "simple.h"

using namespace std;

typedef Eigen::VectorXd vec;

int main() {

cout << "test block in Eigen" << endl;
vec v(3);
v[0] = 0.0;
v[1] = 1.0;
v[2] = 3.0;

cout << "v = " << endl << v << endl;
cout << "v.head(1) = " << endl  << v.head(1) << endl;
cout << "v.tail(1) = " << endl  << v.tail(1) << endl;
cout << "v.segment(1,2) = " << endl  << v.segment(1,2) << endl;

return 0;
}
