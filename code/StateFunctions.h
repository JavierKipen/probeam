#pragma once
#include "Params.h"
#include "State.h"


unsigned int getCountN(State& s);		//Returns count of total N
unsigned int getCountK(State& s);		//Returns count of total K
bool finishedSequencing(State& s);		//Returns if K==0 for every element
bool isEqual(State& s1, State& s2);		//Compares N,K and R to check if they are equal
bool isNEqual(State& s1, unsigned int Np[N_COLORS]);		//Compares N
bool isNEqual(unsigned int N1[N_COLORS], unsigned int N2[N_COLORS]);		//Compares N
float pDyeAttached(State& s, unsigned int i);	//Returns probability of dye attached.
float logpObs(State& s, float obs[N_COLORS]);	//Returns log probability of observing the data given the state
float estExpTerm(unsigned int K[N_COLORS], float obs_norm[N_COLORS]);	//Estimates exponent term
State* decK(State *s, unsigned int i);	//Decreases K at index i(MODIFIES STATE)
State* decN(State* s, unsigned int i);	//Decreases K at index i(MODIFIES STATE)
State* detach(State* s);					//Sets K to 0(MODIFIES STATE)
State* setK(State* s, unsigned int Kp[N_COLORS]);	//Sets K to a given value(MODIFIES STATE)
State* setN(State* s, unsigned int Np[N_COLORS]);	//Sets N to a given value(MODIFIES STATE)
State* setNewDyeSeqsIdxs(State* s, unsigned int * dyeSeqsIdxs, unsigned int newDyeSeqsIdxsCount);	//Copies new dye seqs idxs
State* remChar(State* s, char ch); //Adds to the removed char sequence the given character

State* copyState(State *dest, State &orig);		//Copies orig state to dest

float getStdState(State& S, unsigned int stdIdx);
float getMeanState(State& S, unsigned int meanIdx);
float probObsSingleColor(float obs, unsigned int k);
float logpObs(unsigned int K[N_COLORS], float obs[N_COLORS]);
unsigned int KDist(unsigned int K1[N_COLORS], unsigned int K2[N_COLORS]);
