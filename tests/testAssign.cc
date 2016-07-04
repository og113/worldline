// quick test of assign for assigning arrays to vectors

#include <vector>
#include <string>
#include <iostream>

using namespace std;

int main() {

// hello
cout << "testAssign: " << endl;

// defining
# define length 4
vector<string> v(length);
string a[length] = {"a","b","c","d"};

// assigning
v.assign(a,a+length);

// printing
for (int j=0; j<length; j++) {
	cout << "a[" << j << "] = " << a[j] << ", ";
	cout << "v(" << j << ") = " << v[j] << endl;
}

return 0;
}
