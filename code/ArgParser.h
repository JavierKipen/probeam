#pragma once
#include <iostream>
#include <string>
using namespace std;

class ArgParser
{
public:
	ArgParser(int argc, char* argv[]);
	unsigned int nBeam;
	string folder_path;
	bool parsedCorrectly;

};

