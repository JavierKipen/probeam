#include "ArgParser.h"

#define DEF_NBEAM 7
// #define DEF_PATH "C:/Users/JK-WORK/Desktop/probeam/probeam/data/NormDatasets/1000Prot/"

#define DEF_PATH "C:/Users/JK-WORK/Documents/modifWhatprot/Own/HMM_modif/Datasets/ForPaper/20000Prot/"

ArgParser::ArgParser(int argc, char* argv[])
{
	nBeam = DEF_NBEAM;
	parsedCorrectly = true;
	if (argc == 1) //No argument provided
		folder_path = DEF_PATH;
	else if (argc == 2)
		folder_path = argv[1];
	else if (argc == 4)
	{
		folder_path = argv[1];
		nBeam = stoi(string(argv[3])); //Assumes argv[2] was -b. 
	}
	else
		parsedCorrectly = false; //Could check more failure cases
}
