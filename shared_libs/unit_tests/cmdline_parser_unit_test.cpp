#include <iostream>
#include <string>
#include "cmdline_parser.h"

#define UNSIGNED_INPUT 31
#define DOUBLE_INPUT 23.456
#define CHAR_INPUT 'c'
#define STRING_INPUT "string"
#define N 3

int main(int argc, char* argv[]) {
	std::cout << "Testing cmdlineParser: ";
	
	unsigned int a;
	double b[N];
	char c;
	std::string s;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<unsigned>("a", &a, 1000, "a [unsigned int]");
    parser.addOptArrayParameter<double>("b", b, 23, "b [double []]", N);
    parser.addOptParameter<char>("c", &c, 'c', "c [char]");
    parser.addOptParameter<std::string>("s", &s, "s", "s [std::string]");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return false;
	
	bool result = true;
	if(result) if(a != UNSIGNED_INPUT) result = false;
	if(result) for(unsigned int i = 0; i < N; i++) if(b[i] != DOUBLE_INPUT) result = false;
	if(result) if(c != CHAR_INPUT) result = false;
	if(result) if(s.compare(STRING_INPUT)) result = false;
	
	std::string res = result ? "passed" : "failed";
	std::cout << res << std::endl << std::endl;
}