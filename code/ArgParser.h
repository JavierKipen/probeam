#pragma once
#include <iostream>
#include <string>
#include "Params.h"
using namespace std;

class ArgParser
{
public:
	ArgParser(int argc, char* argv[]);
	unsigned int nBeam,nSparsity;
	string folder_path;
	bool parsedCorrectly;
	float cutoffTh;

};

