#pragma once
#include <iostream>
#include <fstream>
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
	string dyeSeqsPath, scoresProbPath, scoreIdsPath;
	string dataPath;
	string basePath;
	string trueLabelsPath;

	ofstream scoresProbFile, scoresIdsFile;
	bool initOk;

	vector<vector<float>> reads; //features will be flattened to simplify encapsulation.
	vector<unsigned int> trueIDs;

	vector<string> dyeSeqs;
	vector<unsigned int> dyeSeqsIdxs;
	vector<unsigned int> dyeSeqsCounts;
	map<unsigned int, unsigned int> dyeSeqsCountsMap;
	DataIO(string folderPath);
	~DataIO();
	
	void getDyeSeqsInfo(void);
	void loadReads(void);
	void loadReads(unsigned int limit);
	void pushScores(vector<unsigned int> &scoresIdxs, vector<float>& scoresProbs);
	void savePredictions(string Path, vector<unsigned int> yPred, vector<float> yPredProb); //Saves the prediction of a recorder into a csv file
	void createMap();
};

