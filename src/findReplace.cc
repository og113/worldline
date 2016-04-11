/*
	findReplace
		cannot get to grips with grep and emacs so writing this to find and replace with wildcards
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char** argv) {

string fileIn = "temp/fgamma2.txt";
string fileOut = "temp/fgamma3.txt";
string in = "Dot";
string out = "dot";

// getting argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("fileIn")==0) fileIn = (string)argv[2*j+2];
		else if (id.compare("fileOut")==0) fileOut = (string)argv[2*j+2];
		else {
			cerr << "argv id " << id << " not understood" << endl;
			return 1;
		}
	}
}

ifstream is;
is.open(fileIn.c_str());
ofstream os;
os.open(fileOut.c_str());

string line;
stringstream ss;

while (is.getline(line)) {
	ss(line);
}

is.close();
os.close();


return 0;
}
