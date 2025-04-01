#include "ArgParser.h"


// #define DEF_PATH "C:/Users/JK-WORK/Desktop/probeam/probeam/data/NormDatasets/1000Prot/"

//#define DEF_PATH "C:/Users/JK-WORK/Documents/modifWhatprot/Own/HMM_modif/Datasets/ForPaper/20000Prot/"

#define DEF_PATH "C:/Users/JK-WORK/Desktop/DatasetsProtInf/20642_Prot/whatprot/"
#define DEF_SPARSITY 100

ArgParser::ArgParser(int argc, char* argv[])
{
	nBeam = DEF_NBEAM;
	cutoffTh = ZSCOREDEFAULT;
	nSparsity = DEF_SPARSITY;
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
	else if (argc == 6)
	{
		folder_path = argv[1];
		nBeam = stoi(string(argv[3])); //Assumes argv[2] was -b. 
		nSparsity = stoi(string(argv[5])); //Assumes argv[4] was -n. 
	}
	else
		parsedCorrectly = false; //Could check more failure cases
}
