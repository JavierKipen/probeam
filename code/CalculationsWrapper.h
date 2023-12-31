#pragma once
#include <iostream>
#include <string>
#include <map>
#include <array>
#include <vector>
#include <unordered_map>
#include "InfoForEdmanDegradation.h"
#include "InitStates.h"


using namespace std;

vector<unsigned int> argsortf(const vector<float>& vf);

class CalculationsWrapper
{
public:
	CalculationsWrapper();
	void init(vector<string>& dyeSeqs, vector<unsigned int>& dyeSeqsIdx, vector<unsigned int>& relCounts, unsigned int nBeam);
	void getFirstMostLikelyStates(vector<State>* outMostLikely, vector<float>* outProbsNorm, float obs[N_COLORS]);
	//void getObsLogProbs(vector<float>* outLogProbs, vector<StateRed>& auxStates, float obs[N_COLORS]);
	void getInfoForEdman(vector<State>& s, unsigned int nStates);
	pair<unsigned int, float> getMostProbDyeSeqIdx(vector<State>& finalStates, vector<float>& finalStatesLogProbs, unsigned int currNStates);
	vector<array<unsigned int, N_COLORS>> KDyeLoss; //Ks to try in the transition;
	vector<float> KProbsDyeLoss; //Probs of Ks to try in the transition;
	string internalDyeSeqIdToStr(unsigned int id);
	void clear();

	vector<IFED> infosEdman;
	float calcDyeLossProb(unsigned int Kinit[N_COLORS], unsigned int Kend[N_COLORS]);
	void orderCompTimesIFED();
private:
	map<array<unsigned int, N_COLORS>, float> IFEDTimeOfState;
	void initMapTimeIFED();
	array<unsigned int, N_COLORS> getStateOriginalN(State& s);
	

	float n_peptides;
	InitStates is; // To obtain the most likely first states.
	void reformatDyeSeqs(vector<string>& dyeSeqs, vector<unsigned int>& dyeSeqsIdx, vector<unsigned int>& dyeSeqsCounts); //formats to chunk of dye sequences togheter en memory, indexes and chunk of probabilities in memory.
	unsigned int nDyeSeqs;
	float l1normMaxDyeLoss;
	vector<char> dyeSeqsTogheter; //All dye sequences concatenated in memory. To read dye sequence 1, dyeSeqsTogheter[dyeSeqsStartIdxsInMem[i]];
	vector<unsigned long> dyeSeqsStartIdxsInMem; //Saves ints of where dye seqs start
	vector<float> relProbs;
	vector<unsigned int> dyeSeqsIdxOUT; //Dye sequences idxs used in the classification. (in the code the idx represent the sapce on the vector
	vector<unsigned int> dyeSeqsCounts;
	
	vector<unsigned int> dyeSeqsOut; //Vector to calculate the final peptide probs
	vector<float> dyeSeqsOutProb;
	float dyeSeqsProbRelOut[N_MAX_DYESEQS_IN_STATE]; //Variables for decoding most likely output
	unsigned int dyeSeqsProbRelOutCount;
	void getRelProbs(State& s);
	
	//To reorder in memory.
	void reOrderDyeSeqs(vector<string>& dyeSeqs, vector<unsigned int>& peptideIdx, vector<unsigned int>& relCounts);
	array<unsigned int, N_COLORS> getNDyeSeq(string& dyeSeq);

};