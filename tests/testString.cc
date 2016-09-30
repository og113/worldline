#include <string>
#include <iostream>
#include <vector>
#include "simple.h"

using namespace std;

int main() {

string line = " c abc d rentifukous 123 -what is -cat 7m_a7 what!!!";

vector<string> words = splitString(line," ");

cout << line << endl;
cout << "number of words = " << words.size() << endl;
for (uint j=0; j<words.size(); j++) {
	cout << words[j] << endl;
}

return 1;
}
