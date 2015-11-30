# include <iostream>

struct scope {
	enum enums{a,b};
};

int main() {

if( __cplusplus == 201103L ) std::cout << "C++11\n" ;
else if( __cplusplus == 199711L ) std::cout << "C++98\n" ;
else std::cout << "pre-standard C++\n" ;

scope::enums e = scope::b;
std::cout << std:: endl << e << std:: endl;

return 0;
}
