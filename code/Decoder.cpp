#include "Decoder.h"
#include <cmath> 

#define pNotDetach ((float)(1 - D))
#define pNotDetachEdmanFails ((float)(1 - D) * E)
#define pNotDetachEdmanWorks ((float)(1 - D) * (1 - E))

#define ZSCOREDEFAULT 5

Decoder::Decoder()
{
	nBeam = NBEAM_DEFAULT;
}

Decoder::Decoder(unsigned int ExtnBeam)
{
	nBeam = ExtnBeam;
	
}
void Decoder::init(vector<string> dyeSeqs, vector<unsigned int> dyeSeqsIdx, vector<unsigned int> relCounts)
{
	cw.init(dyeSeqs, dyeSeqsIdx, relCounts, nBeam);
	auxStates.reserve(N_CAND_STATE_RESERVE);
	auxStatesProb.reserve(N_CAND_STATE_RESERVE); 

	for (unsigned int i = 0; i < N_FEATURES_PER_COL; i++) {
		vector<State> auxState(nBeam);
		vector<float> auxProb(nBeam);
		mostLikelyStates.push_back(auxState);
		mostLikelyStatesProbNorm.push_back(auxProb);
	}
	obsLogProb.reserve(N_CAND_STATE_RESERVE);
	currT = 0;
	noNextStates = false;
	zScoreTh = ZSCOREDEFAULT;
	PosNextKs.reserve(100); //Ks where next states can finish.
	PosNextKsObsProbLog.reserve(100); //Ks where next states can finish.
	vector<StateRed> auxStates; //States kept for the recursion
	currNStates = nBeam;
	prevStateIndex.reserve(N_CAND_STATE_RESERVE); //Index that says from which state was generated, to recover dyeseqs information
}

pair<unsigned int, float> Decoder::decode(float rad[N_FEATURES_PER_COL][N_COLORS])
{
	cw.getFirstMostLikelyStates(&(mostLikelyStates[0]), &(mostLikelyStatesProbNorm[0]), rad[0]);
	for (currT = 1; currT < N_FEATURES_PER_COL; currT++)
	{
		recursion(rad[currT]);
		if (earlyFinish()) //Checks if all the states are already at K=0. In that case the decoding is already finished.
			break;
	}
	if (currT == N_FEATURES_PER_COL) //If there was no break, takes the last states for picking the last decoding.
		currT = N_FEATURES_PER_COL - 1;
	if (noNextStates) //When no states were longer found, uses the previous best states to predict the peptide.
		currT--;
	pair<unsigned int, float> retVal = getMostProbDyeSeqIdx();
	
	clear();
	return retVal;
}

void Decoder::recursion(float obs[N_COLORS])
{
	calcTransitionProbs(obs); //Obtains next possible states and the probability to transition to them
	calcObservationProbs(obs); //Multiplies the transition probs by the obs probs
	keepBestStates(); //Out of all the aux states keeps the bests.
	auxStates.clear(); auxStatesProb.clear(); obsLogProb.clear(); prevStateIndex.clear(); //Clears auxiliar states for next recursion.
}

void Decoder::calcTransitionProbs(float obs[N_COLORS])
{
	unsigned int K_aux[N_COLORS];
	cw.getInfoForEdman(mostLikelyStates[currT - 1], currNStates);
	getBestKObs(obs); //Gets the most likely K, which is used for selecting which possible outputs can we have.
	getPosNextK(obs); //Then keeps the Ks that are close to the observation
	for (unsigned int i = 0; i < currNStates; i++)
	{
		State& Si = mostLikelyStates[currT - 1][i];
		IFED& info = cw.infosEdman[i];
		
		float stateLogProb = mostLikelyStatesProbNorm[currT - 1][i];
		float stateProb = (float)exp(stateLogProb);
		for (unsigned int Kidx=0; Kidx < PosNextKs.size(); Kidx++)
		{
			if (isK1leqK2elem(PosNextKs[Kidx].data(),Si.K))
			{
				calcAndPushPosStates(Si, Kidx, stateProb, info,i);
			}
		}
	}
	
}

void Decoder::calcAndPushPosStates(State& Si, unsigned int posNextKIdx, float stateProb, IFED& infoEd,unsigned int idxPrev)
{
	float pDetach = D;
	float pNoDetach = (1 - D);
	float pSuccessRem = pNoDetach * (1 - E); //Did not detach and Edman process not failed.
	array<unsigned int, N_COLORS> posNextK = PosNextKs[posNextKIdx];
	float dyeLossProb= cw.calcDyeLossProb(Si.K, posNextK.data());
	array<unsigned int, N_COLORS> auxK= posNextK;
	StateRed auxState;
	setK(copyState(&auxState, Si), posNextK.data()); //Copies the state and sets the K to the value that is being analyzed
	
    float noRemProb = stateProb * (1 - D) * E * (dyeLossProb); // Did not detach, fails edman cycle * dye loss prob
	
	if (isNZero(Si.K))
		appendCandidateState(&auxState, stateProb, idxPrev, posNextKIdx); //When K=0, state goes to self. Any other effect does not matter.
	else
	{
		//No removal
		if (isNZero(posNextK.data())) //If it was 0, we add the detachment prob.
			noRemProb += stateProb * D;
		appendCandidateState(&auxState, noRemProb, idxPrev, posNextKIdx); //State due only to dye loss, Edman cycle fails
		//Removal of a "."
		if (infoEd.nonLumCanBeRemoved)
		{
			appendCandidateState(remChar(&auxState, '.'), stateProb * (1 - D) * (1 - E) * dyeLossProb * infoEd.pRemNonLum, idxPrev, posNextKIdx); //State due only to dye loss, Edman cycle fails
			omitRemoval(&auxState); //Resets aux var
		}
		//Removal of a dye
		for (unsigned int i = 0; i < N_COLORS; i++)
		{
			if (infoEd.dyeCanBeRemoved[i])
			{
				float pDA = pDyeAttached(Si, i);
				if (pDA != 0 && Si.K[i] > posNextK[i]) //Considering the case that the dye was attached
				{  //And posNextK has less K than the original state
					auxK[i] += 1;
					if (isK1leqK2elem(Si.K, auxK.data()))
					{
						appendCandidateState(remChar(decN(&auxState, i), i + '0'), stateProb * (1 - D) * (1 - E) * cw.calcDyeLossProb(Si.K, auxK.data()) * pDA * infoEd.pRemDye[i], idxPrev, posNextKIdx);
						omitRemoval(incN(&auxState, i));
					}
					//Decrease S.Ni--,S.R+=i
				}
				if (pDA != 1) //Dye not attached
				{
					appendCandidateState(remChar(decN(&auxState, i), i + '0'), stateProb * (1 - D) * (1 - E) * dyeLossProb * (1 - pDA) * infoEd.pRemDye[i], idxPrev, posNextKIdx);
					omitRemoval(incN(&auxState, i));
				}
			}
		}
	}
}

bool Decoder::earlyFinish() //Finish sequencing if it does not have next possible states or all states are with K=0.
{
	bool allZeros = true;
	for (unsigned int i = 0; i < mostLikelyStates[currT].size(); i++)
	{
		allZeros &= isNZero(mostLikelyStates[currT][i].K);
		if (!allZeros)
			break;
	}
	return allZeros||noNextStates;
}

void Decoder::calcObservationProbs(float obs[N_COLORS])
{
	for (unsigned int i = 0; i < auxStates.size(); i++)
		auxStatesProb[i] = log(auxStatesProb[i]) + obsLogProb[i]; // Log prob of the state + log prob of observation.
}
void Decoder::keepBestStates()
{
	
	unsigned int currStateIdx;
	unsigned int stateCount = auxStates.size();
	unsigned int statesToUse = nBeam;
	if (stateCount == 0)
		noNextStates = true;
	else
	{
		vector<unsigned int> idxsSorted = argsortf(auxStatesProb);
		if (stateCount < nBeam) // If we get less states in one t, we only use those.
			statesToUse = stateCount;
		for (unsigned int i = 0; i < statesToUse; i++)
		{
			currStateIdx = idxsSorted[i];
			copyState(&mostLikelyStates[currT][i], auxStates[currStateIdx]);
			retrieveDyeSequences(&mostLikelyStates[currT][i], currStateIdx);
			if (i == 0)
				mostLikelyStatesProbNorm[currT][i] = auxStatesProb[currStateIdx];
			else
				mostLikelyStatesProbNorm[currT][i] = auxStatesProb[currStateIdx] - mostLikelyStatesProbNorm[currT][0];
		}
		mostLikelyStatesProbNorm[currT][0] = 0; //Normalizes first value at the very end
		currNStates = statesToUse;
	}
	
}

pair<unsigned int, float> Decoder::getMostProbDyeSeqIdx()
{
	return cw.getMostProbDyeSeqIdx(mostLikelyStates[currT], mostLikelyStatesProbNorm[currT]);
}

void Decoder::appendCandidateState(StateRed* s, float prob, unsigned int idxPrev, unsigned int posNextKIndex)
{
	bool found = false; //Pushes state if it was not here before, if it was adds the probability to transition to that state.
	if (prob < 0)
		cout << "This shouldnt happen";
	for (unsigned int i = 0; i < auxStates.size(); i++) 
	{
		if (isEqual(auxStates[i], *s))
		{
			auxStatesProb[i] += prob; //If state was before, adds probabilities
			found = true;
			break;
		}
	}
	if (found == false)
	{
		auxStates.push_back(*s);
		auxStatesProb.push_back(prob);
		prevStateIndex.push_back(idxPrev);
		obsLogProb.push_back(PosNextKsObsProbLog[posNextKIndex]);
	}
}

void Decoder::getBestKObs(float obs[N_COLORS])
{
	for (unsigned int i = 0; i < N_COLORS; i++)
		currBestK[i] = getBestKi(obs[i]);
}

void Decoder::getPosNextK(float obs[N_COLORS])
{
	PosNextKs.clear();
	PosNextKsObsProbLog.clear();
	PosNextKs.push_back(currBestK); //Pushes the best K
	PosNextKsObsProbLog.push_back(logpObs(currBestK.data(), obs));
	//Explore recursively neighbours
	array<unsigned int, N_COLORS> aux;
	aux = currBestK;
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		if (aux[i] > 0)
		{
			aux[i] -= 1;
			recursivePosNextKs(obs, aux);
			aux[i] += 1;
		}
		aux[i] += 1;
		recursivePosNextKs(obs, aux);
		aux[i] -= 1;
	}

}

void Decoder::recursivePosNextKs(float obs[N_COLORS],array<unsigned int, N_COLORS> & origK)
{
	if (approxDistZScore(origK.data(), obs) < zScoreTh)
	{
		bool found = false;
		for (unsigned int i = 0; i < PosNextKs.size(); i++)
		{
			if (isNEqual(PosNextKs[i].data(), origK.data()))
			{
				found = true;
				break;
			}
		}
		if (found == false)
		{
			PosNextKs.push_back(origK);
			PosNextKsObsProbLog.push_back(logpObs(origK.data(), obs));
			array<unsigned int, N_COLORS> aux;
			aux = origK;
			for (unsigned int i = 0; i < N_COLORS; i++)
			{
				if (aux[i] > 0)
				{
					aux[i] -= 1;
					recursivePosNextKs(obs, aux);
					aux[i] += 1;
				}
				aux[i] += 1;
				recursivePosNextKs(obs, aux);
				aux[i] -= 1;
			}
		}
	}
}

void Decoder::retrieveDyeSequences(State* mostlikelyState, unsigned int selectedAuxStateIdx)
{
	StateRed& s = auxStates[selectedAuxStateIdx];
	unsigned int prevStateIdx = prevStateIndex[selectedAuxStateIdx];
	State& prevState = mostLikelyStates[currT-1][prevStateIdx];
	if (prevState.RCharCount == s.RCharCount) //No removal
	{
		for (unsigned int i = 0; i < prevState.dyeSeqsIdxsCount; i++)
			mostlikelyState->dyeSeqsIdxs[i] = prevState.dyeSeqsIdxs[i];
		mostlikelyState->dyeSeqsIdxsCount = prevState.dyeSeqsIdxsCount;
	}
	else
	{
		char removedaa = s.R[s.RCharCount - 1];
		if (removedaa == '.')
		{
			unsigned int totalCount = cw.infosEdman[prevStateIdx].dyeSeqsIdxsCountDot;
			unsigned int* dyeSeqs = cw.infosEdman[prevStateIdx].dyeSeqsIdxsDot;
			for (unsigned int i = 0; i < totalCount; i++)
				mostlikelyState->dyeSeqsIdxs[i] = dyeSeqs[i];
			mostlikelyState->dyeSeqsIdxsCount = totalCount;
		}
		else
		{
			unsigned int dyeRem = removedaa - '0';
			unsigned int totalCount = cw.infosEdman[prevStateIdx].dyeSeqsIdxsCount[dyeRem];
			unsigned int* dyeSeqs = cw.infosEdman[prevStateIdx].dyeSeqsIdxs[dyeRem];
			for (unsigned int i = 0; i < totalCount; i++)
				mostlikelyState->dyeSeqsIdxs[i] = dyeSeqs[i];
			mostlikelyState->dyeSeqsIdxsCount = totalCount;
		}
	}
}

void Decoder::clear()
{
	auxStates.clear();
	auxStatesProb.clear();
	cw.clear();
	currT = 0;
	noNextStates = false;
	currNStates = nBeam;
	/*for (unsigned int i = 0; i < N_FEATURES_PER_COL; i++) {
		mostLikelyStates[i].clear();
		mostLikelyStatesProbNorm[i].clear();
	}*/
}

//Useful functions outside of class
