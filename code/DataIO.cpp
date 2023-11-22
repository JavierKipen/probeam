#include "DataIO.h"
#include <fstream>
#include <sstream>


using namespace std;


DataIO::DataIO(string folderPath)
{
	dyeSeqsPath = folderPath + "dye-seqs.tsv";
	dataPath = folderPath + "radiometries.tsv";
	trueLabelsPath = folderPath + "true-ids.tsv";
	basePath = folderPath ;
	//trueIDs.reserve(10000); //Usually datasets will have around 10k reads
	//trueIDs.reserve(10000); //Usually datasets will have around 10k reads
	getDyeSeqsInfo();
	initOk = true; //Should check that everything exists.
}

void DataIO::getDyeSeqsInfo(void)
{
	fstream dyeSeqsFile(dyeSeqsPath, ios::in);
	string auxLine;
	unsigned int dyeSeqsCount;
	if (!dyeSeqsFile.is_open())
		cout << "dyeSeqs file not found!";
	else 
	{
		getline(dyeSeqsFile, auxLine); //reads first line, says number of colors. We dont use it.
		getline(dyeSeqsFile, auxLine); //reads second line, says how many dye seqs we have. We use it to reserve memory in our vectors.
		dyeSeqsCount = stoi(auxLine);
		dyeSeqs.reserve(dyeSeqsCount); dyeSeqsIdxs.reserve(dyeSeqsCount); dyeSeqsCounts.reserve(dyeSeqsCount); //Reserving memory

		while (getline(dyeSeqsFile, auxLine)) //Now starts parsing data from the dye seq file.
		{
			stringstream ss(auxLine);
			vector<string> substrings;
			while (ss.good()) { //Gets the 3 substrings for each line
				string substr;
				getline(ss, substr, '\t');
				substrings.push_back(substr);
			}
			dyeSeqs.push_back(substrings[0]);
			dyeSeqsCounts.push_back(stoi(substrings[1]));
			dyeSeqsIdxs.push_back(stoi(substrings[2]));
		}
		dyeSeqsFile.close();

	}
}

void DataIO::loadReads(void)
{
	fstream dataFile(dataPath, ios::in);
	fstream trueIDsFile(trueLabelsPath, ios::in);
	vector<float> auxVect(N_FEATURES_TOTAL);
	string auxLine, feature;
	unsigned int featureCount = 0;
	unsigned int readCount = 0;
	if (!dataFile.is_open())
		cout << "Radiometries file not found!";
	else
	{
		getline(dataFile, auxLine); //reads first line, N Edman+1
		getline(dataFile, auxLine); //reads Second line, N colors
		getline(dataFile, auxLine); //reads third line, number of reads. This one is useful:
		reads.reserve(stoi(auxLine));
		for (unsigned int i = 0; i < stoi(auxLine); i++)
			reads.push_back(auxVect);
		trueIDs.resize(stoi(auxLine));
		while (getline(dataFile, auxLine))
		{
			unsigned int index;
			featureCount = 0;
			stringstream ss(auxLine);
			while (ss.good()) {
				getline(ss, feature, '\t');
				//index = (featureCount % N_COLORS) * N_FEATURES_PER_COL + featureCount / N_COLORS;
				reads[readCount][featureCount] = stof(feature);
				featureCount++;
			}
			readCount++;
		}
		dataFile.close();
		if (!trueIDsFile.is_open())
			cout << "True IDs file not found!";
		else
		{
			readCount = 0;
			getline(trueIDsFile, auxLine); //reads first line, number of ids.
			while (getline(trueIDsFile, auxLine))
				trueIDs[readCount++] = stoi(auxLine);
			trueIDsFile.close();
		}
	}
}

void DataIO::loadReads(unsigned int limit)
{
	fstream dataFile(dataPath, ios::in);
	fstream trueIDsFile(trueLabelsPath, ios::in);
	string auxLine, feature;
	vector<float> auxVect(N_FEATURES_TOTAL);
	unsigned int featureCount = 0;
	unsigned int readCount = 0;
	if (!dataFile.is_open())
		cout << "Radiometries file not found!";
	else
	{
		getline(dataFile, auxLine); //reads first line, N Edman+1
		getline(dataFile, auxLine); //reads Second line, N colors
		getline(dataFile, auxLine); //reads third line, number of reads. This one is useful:
		reads.reserve(limit);
		for (unsigned int i = 0; i < limit; i++)
			reads.push_back(auxVect);
		trueIDs.resize(limit);
		
		while (readCount<limit)
		{
			unsigned int index;
			featureCount = 0;
			getline(dataFile, auxLine);
			stringstream ss(auxLine);
			while (ss.good()) {
				getline(ss, feature, '\t');
				//index = (featureCount % N_COLORS) * N_FEATURES_PER_COL + featureCount / N_COLORS;
				reads[readCount][featureCount] = stof(feature);
				featureCount++;
			}
			readCount++;
		}
		dataFile.close();
		if (!trueIDsFile.is_open())
			cout << "True IDs file not found!";
		else
		{
			readCount = 0;
			getline(trueIDsFile, auxLine); //reads first line, number of ids.
			while (readCount < limit && getline(trueIDsFile, auxLine))
				trueIDs[readCount++] = stoi(auxLine);
			trueIDsFile.close();
		}
	}
}

void DataIO::savePredictions(string path, vector<unsigned int> yPred, vector<float> yPredProb)
{
	ofstream myfile;
	myfile.open(path);
	myfile << "radmat_iz, best_pep_iz, best_pep_score\n";
	for (unsigned int i = 0; i < yPred.size(); i++)
	{
		myfile << to_string(i) << "," << to_string(yPred[i]) << "," << to_string(yPredProb[i]) << "\n" ;
	}
	myfile.close();
}

void DataIO::createMap()
{
	for (unsigned int i = 0; i < dyeSeqs.size(); i++)
		dyeSeqsCountsMap.insert(pair<unsigned int, unsigned int>(dyeSeqsIdxs[i], dyeSeqsCounts[i]));
}

