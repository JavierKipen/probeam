#pragma once
#include <iostream>
#include "State.h"
#include "CalculationsWrapper.h"
#include "StateFunctions.h"
#include <utility> 

#define NBEAM_DEFAULT 10
#define N_CAND_STATE_RESERVE 1000 //Number of states reserved for candidate states finding


class Decoder
{

public:
	Decoder();
	Decoder(unsigned int nBeam);
	unsigned int nBeam;
	void init(vector<string> dyeSeqs, vector<unsigned int> dyeSeqsIdx, vector<unsigned int> relCounts); //Initializes internal variables
	pair<unsigned int, float> decode(float rad[N_FEATURES_PER_COL][N_COLORS]); //Decodes a given read. Returns id of dye sequence and probability
	CalculationsWrapper cw;
private:
	
	vector<vector<State>> mostLikelyStates; //States kept for every time T
	vector<vector<float>> mostLikelyStatesProbNorm; //Prob of states kept for every time T
	vector<State> auxStates; //States kept for the recursion
	vector<float> auxStatesProb; //Prob of states kept for the recursion
	vector<float> auxObsLogProb; //Used for the observation probabilities in the aux calc
	void recursion(float obs[N_COLORS], unsigned int t);
	void calcTransitionProbs(float obs[N_COLORS], unsigned int  t);
	void transProbSuccesfulRemoval(State &S, float initProb); //InitProb contains already the probability of dye loss + prob not detached.
	void calcObservationProbs(float obs[N_COLORS]);
	void keepBestStates(unsigned int t);
	pair<unsigned int, float> getMostProbDyeSeqIdx(); //This function gets the most likely states at the last time and outputs the most likely dye sequences.
	void clear(); //To reset parameters after each decoding.
	void appendCandidateState(State& s, float prob);

};

