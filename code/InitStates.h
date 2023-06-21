#pragma once
#include <iostream>
#include <string>
#include <map>
#include <array>
#include <vector>
#include "InfoForEdmanDegradation.h"
#include "StateFunctions.h"


using namespace std;


vector<unsigned int> argsortf(const vector<float>& vf);
unsigned int nChoosek(unsigned int n, unsigned int k);
void highestKValsInArray(float* array, unsigned long n, unsigned int k, unsigned long* outBestIdx);

typedef struct {
	unsigned int stateIndex, K;
	float logProb;
} InitStateInfo; //For ordering in smart search


class InitStates
{
public:
	InitStates();
	void init(vector<char> *dyeSeqsTogheter, vector<unsigned long> *dyeSeqsStartIdxsInMem, vector<float> *relProbs,unsigned int nBeamIn);
	void getFirstMostLikelyStates(vector<State>* outMostLikely, vector<float>* outProbsNorm,  float obs[N_COLORS]);
	void clear();
private:
	//Information of the dye sequences
	unsigned int nDyeSeqs;
	vector<char> *dyeSeqsTogheter; //All dye sequences concatenated in memory. To read dye sequence 1, dyeSeqsTogheter[dyeSeqsStartIdxsInMem[i]];
	vector<unsigned long> *dyeSeqsStartIdxsInMem; //Saves ints of where dye seqs start
	vector<float> *relProbs;
	vector<unsigned int> dyeSeqsIdxOUT; //Dye sequences idxs used in the classification. (in the code the idx represent the sapce on the vector
	//Ideal init states and their probs
	vector<float> initIdealStatesProb; 
	vector<float> initIdealStatesLogProb;
	vector<State> initIdealStates;
	void createInitIdealStates();
	void appendForIdInSt(unsigned int currN[N_COLORS], unsigned int dyeSeqIdx, float prob);

	//For Hardcore init states (checking every single probability with every posible state.
	vector<float> initMeansH, initStdsH, initLogProbTermH, initStatesLogProbH, initStatesLogProbHCopy; //Variables for hardcore init most likely states.
	vector<State> initStates;
	void extendToInitStates();
	void hardcoreMostLikelyInitStates(vector<State>* outMostLikely, vector<float>* outProbsNorm, unsigned int nBeam, float obs[N_COLORS]);
	void precomputeForHardcoreInitStates();

	//For smart init states
	unsigned int nBeam;
	unsigned int maxNis[N_COLORS];
	unsigned int shiftsMin[N_COLORS];
	unsigned int KBest[N_COLORS];
	vector<float> maxIncDyeMiss;
	vector<map<unsigned int, float>> PdyeMissInit; //For each ideal states has the dye miss probability of each state.
	map<unsigned int, float> obsProbMap; //Map that contains the log of gaussian evaluated, so it doesnt have to be calculated again
	vector<unsigned int> bestStatesIdx;
	vector<InitStateInfo> bestStatesInfo;
	void getMaxNi();
	void smartSearchStates(float obs[N_COLORS]);
	unsigned int NorKToUnsInt(unsigned int NorK[N_COLORS]);
	void KuintToVector(unsigned int NorK[N_COLORS], unsigned int baseK);
	void generatePosKsInit();
	vector<unsigned int> orderByClosestStates(float obs[N_COLORS]);
	void getBestK(float obs[N_COLORS]);
	void smartMostLikelyInitStates(vector<State>* outMostLikely, vector<float>* outProbsNorm,  float obs[N_COLORS]);
	void preLoadBestStates(float obs[N_COLORS]);
	float getObsLogProb(float obs[N_COLORS], unsigned int K_uint);
	void getBestKState(unsigned int K[N_COLORS], unsigned int N[N_COLORS]);
	void pushInitStateInfo(InitStateInfo &aux);
	void exploreOtherKs(InitStateInfo &aux, unsigned int sKBest[N_COLORS], float obs[N_COLORS]);
	void formatOutput(vector<State>* outMostLikely, vector<float>* outProbsNorm);
	
};

