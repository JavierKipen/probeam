#include <iostream>
#include <string>
#include <chrono>
#include "State.h"
#include "DataIO.h"
#include "Decoder.h"
#include "ArgParser.h"


using namespace std;

void checkOnlyOneRead(Decoder &dec);
void checkReducedData(Decoder& dec, DataIO &dataIO,unsigned int red);
void checkWholeDataset(Decoder& dec, DataIO &dataIO);
void timeRunAndSave(Decoder& dec, DataIO& dataIO);

int main(int argc, char* argv[])
{
	ArgParser argParser(argc, argv);
	if (argParser.parsedCorrectly)
	{
		DataIO dataIO(argParser.folder_path);
		if (dataIO.initOk)
		{
			Decoder decoder(argParser.nBeam);
			decoder.init(dataIO.dyeSeqs, dataIO.dyeSeqsIdxs, dataIO.dyeSeqsCounts);
			timeRunAndSave(decoder, dataIO);
			//checkReducedData(decoder, dataIO, 2000);
			//checkWholeDataset(decoder, dataIO);
		}
		else
			cout << "Error generating model";
	}
	else
		cout << "Error in argument parsing";
	return 0;
}



void checkOnlyOneRead(Decoder& dec)
{	 
	float rad[N_FEATURES_PER_COL][N_COLORS];
	//rad[0][0] = 28570; rad[0][1] = 12783; rad[0][2] = -62;
	//rad[1][0] = 27569; rad[1][1] = 10870; rad[1][2] = 27;
	//rad[2][0] = 22102; rad[2][1] = 11733; rad[2][2] = 0.14;
	//rad[3][0] = 19565; rad[3][1] = 9937;  rad[3][2] = 90.25;
	//rad[4][0] = 9848;  rad[4][1] = 11046; rad[4][2] = -50;
	//rad[5][0] = 9166;  rad[5][1] = 10554; rad[5][2] = 15.02;
	//rad[6][0] = 12433; rad[6][1] = 8536; rad[6][2] = 42.6;
	//rad[7][0] = -69.37;rad[7][1] = 12303; rad[7][2] = 9.3;
	//rad[8][0] = 5;     rad[8][1] = -7.3;  rad[8][2] = -3.72;
	//rad[9][0] = 46.07; rad[9][1] = 8.9; rad[9][2] = -7.4;

	 rad[0][0] = 10039; rad[0][1] = 0; rad[0][2] = -62;
	 rad[1][0] = 9600; rad[1][1] = 0; rad[1][2] = 0;
	 rad[2][0] = 0; rad[2][1] = 0; rad[2][2] = 0.14;
	 rad[3][0] = 0; rad[3][1] = 0;  rad[3][2] = 90.25;
	 rad[4][0] = 0;  rad[4][1] = 0; rad[4][2] = -50;
	 rad[5][0] = 0;  rad[5][1] = 0; rad[5][2] = 15.02;
	 rad[6][0] = 0; rad[6][1] = 0; rad[6][2] = 42.6;
	 rad[7][0] = -0; rad[7][1] = 0; rad[7][2] = 9.3;
	 rad[8][0] = 5;     rad[8][1] = -7.3;  rad[8][2] = -3.72;
	 rad[9][0] = 46.07; rad[9][1] = 8.9; rad[9][2] = -7.4;
	 auto start = chrono::high_resolution_clock::now();
	 pair<unsigned int, float> auxRes;
	 auxRes = dec.decode(rad);
	 auto stop = chrono::high_resolution_clock::now();
	 auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	 cout << "Time Decoding one read: "
		 << duration.count() << " microseconds" << endl; //https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
	 cout << to_string(auxRes.first);
}
void checkReducedData(Decoder& dec, DataIO &dataIO,unsigned int red)
{
	unsigned int reduced = red;
	dataIO.loadReads(reduced);
	vector<unsigned int> yPred;
	vector<float> yPredProb;
	yPred.reserve(dataIO.reads.size());
	yPredProb.reserve(dataIO.reads.size());
	pair<unsigned int, float> auxRes;
	auto start = chrono::high_resolution_clock::now();
	for (unsigned int i = 0; i < reduced; i++)
	{
		auxRes = dec.decode((float(*)[3])dataIO.reads[i].data());
		//auxRes = decoder.decode((float(*)[3])dataIO.reads[1].data());
		yPred.push_back(auxRes.first);
		yPredProb.push_back(auxRes.second);
	}
	auto stop = chrono::high_resolution_clock::now();
	unsigned int correctCount = 0;
	for (unsigned int i = 0; i < reduced; i++)
		if (yPred[i] == dataIO.trueIDs[i])
			correctCount++;
	cout << "Accuracy: " + to_string(((float)correctCount) / ((float)reduced)) << endl;
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	cout << "Time Decoding one read: "
		<< duration.count() / reduced << " microseconds" << endl;
	//dataIO.savePredictions(folder_path + "BeamSearchPred10.csv", yPred, yPredProb);
	
}
void checkWholeDataset(Decoder& dec, DataIO& dataIO )
{
	dataIO.loadReads();
	vector<unsigned int> yPred;
	vector<float> yPredProb;
	yPred.reserve(dataIO.reads.size());
	yPredProb.reserve(dataIO.reads.size());
	pair<unsigned int, float> auxRes;
	auto start = chrono::high_resolution_clock::now();
	for (unsigned int i = 0; i < dataIO.reads.size(); i++)
	{
		auxRes = dec.decode((float(*)[3])dataIO.reads[i].data());
		yPred.push_back(auxRes.first);
		yPredProb.push_back(auxRes.second);
	}
	auto stop = chrono::high_resolution_clock::now();
	unsigned int correctCount = 0;
	for (unsigned int i = 0; i < dataIO.reads.size(); i++)
		if (yPred[i] == dataIO.trueIDs[i])
			correctCount++;
	cout << "Accuracy: " + to_string(((float)correctCount) / ((float)dataIO.reads.size())) << endl;
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	cout << "Time Decoding one read: "
		<< duration.count() / dataIO.reads.size() << " microseconds" << endl;
	//dataIO.savePredictions(folder_path + "BeamSearchPred.csv", yPred, yPredProb);
}

void timeRunAndSave(Decoder& dec, DataIO& dataIO)
{
	dataIO.loadReads();
	vector<unsigned int> yPred;
	vector<float> yPredProb;
	yPred.reserve(dataIO.reads.size());
	yPredProb.reserve(dataIO.reads.size());
	pair<unsigned int, float> auxRes;
	auto start = chrono::high_resolution_clock::now();
	for (unsigned int i = 0; i < dataIO.reads.size(); i++)
	{
		//if (i % 10 == 0)
		//	cout << i << endl; //1870 explodes
		auxRes = dec.decode((float(*)[3])dataIO.reads[i].data());
		yPred.push_back(auxRes.first);
		yPredProb.push_back(auxRes.second);
	}
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	cout << "Time Decoding one read: " << (float)duration.count() / (float)dataIO.reads.size() << " microseconds" << endl;
	dataIO.savePredictions(dataIO.basePath + "BeamSearchPred"+ to_string(dec.nBeam)+".csv", yPred, yPredProb);
}