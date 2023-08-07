#pragma once
#include <iostream>
#include "State.h"
#include "CalculationsWrapper.h"
#include "StateFunctions.h"
#include <utility> 

#define NBEAM_DEFAULT 15
#define N_CAND_STATE_RESERVE 1000 //Number of states reserved for candidate states finding



class Decoder
{

public:
	Decoder();
	Decoder(unsigned int nBeam);
	Decoder(unsigned int nBeam, float cutoffTh);
	
	unsigned int nBeam;
	void init(vector<string> dyeSeqs, vector<unsigned int> dyeSeqsIdx, vector<unsigned int> relCounts); //Initializes internal variables
	pair<unsigned int, float> decode(float rad[N_FEATURES_PER_COL][N_COLORS]); //Decodes a given read. Returns id of dye sequence and probability
	CalculationsWrapper cw;
private:
	unsigned int currT,currNStates;
	float zScoreTh;
	bool noNextStates;
	array<unsigned int, N_COLORS> currBestK; //The most likely K for a Xt, used in the transitions.
	vector<array<unsigned int, N_COLORS>> PosNextKs; //Ks where next states can finish.
	vector<float> PosNextKsObsProbLog; //Ks where next states can finish.
	bool earlyFinish();
	void getPosNextK(float obs[N_COLORS]);
	void getBestKObs(float obs[N_COLORS]);
	void calcAndPushPosStates(unsigned int posNextKIdx, unsigned int idxStatePrev);
	vector<StateRed> auxStates; //States kept for the recursion
	vector<unsigned int> prevStateIndex; //Index that says from which state was generated, to recover dyeseqs information
	vector<float> auxStatesProb; //Prob of states kept for the recursion
	vector<float> obsLogProb; //Prob of states kept for the recursion
	void appendCandidateState(StateRed *s, float prob, unsigned int idxPrev, unsigned int posNextKIndex);
	void retrieveDyeSequences(State * mostlikelyState,unsigned int selectedAuxStateIdx);
	void recursivePosNextKs(float obs[N_COLORS], array<unsigned int, N_COLORS>& origK);
	void uponSuccessRem(unsigned int idxStatePrev, unsigned int posNextKIndex, unsigned int dyeIdx);

	vector<vector<State>> mostLikelyStates; //States kept for every time T
	vector<vector<float>> mostLikelyStatesProbNorm; //Prob of states kept for every time T
	

	void recursion(float obs[N_COLORS]);
	void calcTransitionProbs(float obs[N_COLORS]);
	void transProbSuccesfulRemoval(State &S, float initProb); //InitProb contains already the probability of dye loss + prob not detached.
	void calcObservationProbs(float obs[N_COLORS]);
	void keepBestStates();
	pair<unsigned int, float> getMostProbDyeSeqIdx(); //This function gets the most likely states at the last time and outputs the most likely dye sequences.
	void clear(); //To reset parameters after each decoding.
	

};

