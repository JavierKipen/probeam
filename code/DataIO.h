#pragma once
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <map>
#include <array>
#include "Params.h"

using namespace std;

class DataIO
{

public:
	string dyeSeqsPath;
	string dataPath;
	string basePath;
	string trueLabelsPath;
	bool initOk;

	vector<vector<float>> reads; //features will be flattened to simplify encapsulation.
	vector<unsigned int> trueIDs;

	vector<string> dyeSeqs;
	vector<unsigned int> dyeSeqsIdxs;
	vector<unsigned int> dyeSeqsCounts;

	DataIO(string folderPath);
	
	void getDyeSeqsInfo(void);
	void loadReads(void);
	void loadReads(unsigned int limit);

	void savePredictions(string Path, vector<unsigned int> yPred, vector<float> yPredProb); //Saves the prediction of a recorder into a csv file

};

