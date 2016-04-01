#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include"print.h"
#include"simple.h"

using namespace std;

int main() {

cout << "vector assign test: ok if following numbers agree pairwise" << endl;

int array[5]={1,2,3,2,1};
cout << array[2] << endl;

vector<int> v(5);
v.assign(array, array + 5);
cout << v[2] << endl;

number arrayNum[5]={1.0,2.0,3.0,4.0,5.0};
cout << arrayNum[2] << endl;

vector<number> vNum(5);
vNum.assign(arrayNum, arrayNum + 5);
cout << vNum[2] << endl;

string file = "tests/v.dat";
saveVectorCsvAppend(file,vNum);
vNum[0] += 9.0;
saveVectorCsvAppend(file,vNum);
cout << "saved vector to " << file << endl;

cout << "number to string test of precision" << endl;
cout << nts(PI) << endl;
cout << nts(PI,1) << endl;
cout << nts(PI,16) << endl;
cout << nts(PI,20) << endl;
cout << nts(PI,30) << endl;

return 0;
}
