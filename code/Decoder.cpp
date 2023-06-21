#include "Decoder.h"
#include <cmath> 

#define pNotDetach ((float)(1 - D))
#define pNotDetachEdmanFails ((float)(1 - D) * E)
#define pNotDetachEdmanWorks ((float)(1 - D) * (1 - E))



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
	auxObsLogProb.reserve(N_CAND_STATE_RESERVE);
}

pair<unsigned int, float> Decoder::decode(float rad[N_FEATURES_PER_COL][N_COLORS])
{
	cw.getFirstMostLikelyStates(&(mostLikelyStates[0]), &(mostLikelyStatesProbNorm[0]), rad[0]);
	for (unsigned int t = 1; t < N_FEATURES_PER_COL; t++)
		recursion(rad[t], t);
	pair<unsigned int, float> retVal = getMostProbDyeSeqIdx();
	
	clear();
	return retVal;
}

void Decoder::recursion(float obs[N_COLORS],unsigned int t)
{
	calcTransitionProbs(obs,t);
	//vector<float> auxDebugProbTrans;
	//for (unsigned int i = 0; i < auxStatesProb.size(); i++)
		//auxDebugProbTrans.push_back(log(auxStatesProb[i]));
	calcObservationProbs(obs);
	keepBestStates(t);
	auxStates.clear(); auxStatesProb.clear(); auxObsLogProb.clear(); //Clears auxiliar states for next recursion.
}

void Decoder::calcTransitionProbs(float obs[N_COLORS], unsigned int  t)
{
	unsigned int K_aux[N_COLORS];
	for (unsigned int i = 0; i < nBeam; i++)
	{
		State Si = mostLikelyStates[t - 1][i];
		float stateLogProb = mostLikelyStatesProbNorm[t - 1][i];
		float stateProb = (float)exp(stateLogProb);
		if (!finishedSequencing(Si))
		{
			K_aux[0] = Si.K[0]; K_aux[1] = Si.K[1]; K_aux[2] = Si.K[2];
			detach(&Si); //Sets Si.K to zero
			appendCandidateState(Si, stateProb * (float) D ); // Pushes Detached state
			setK(&Si, K_aux); //Restores values.
			cw.obtainKDyeLoss(Si, obs);			//Calculator obtains the possible dye losses and their probs (in cw.KDyeLoss, cw.KProb...)
			for (unsigned int j = 0; j < cw.KDyeLoss.size(); j++)
			{
				setK(&Si, cw.KDyeLoss[j].data()); //Restores values.
				appendCandidateState(Si, stateProb * pNotDetachEdmanFails * cw.KProbsDyeLoss[j]); //Not detached, edman cycle fails;
				transProbSuccesfulRemoval(Si, stateProb * pNotDetachEdmanWorks * cw.KProbsDyeLoss[j]); //Not detached edman cycle worked.
			}
			cw.KDyeLoss.clear(); cw.KProbsDyeLoss.clear();
		}
		else
			appendCandidateState(Si, stateProb);//Probability 1 that the state remains the same!
	}
}

void Decoder::transProbSuccesfulRemoval(State &s, float initProb)
{
	State sn(s); //next state. before modification is always a hard copy of original s.
	cw.getInfoForEdman(s);
	IFED &info = cw.infoEdman;
	float pDA;
	if (info.nonLumCanBeRemoved)
	{
		appendCandidateState(*(setNewDyeSeqsIdxs(remChar(&sn, '.'), info.dyeSeqsIdxsDot,info.dyeSeqsIdxsCountDot)), 
			info.pRemNonLum * initProb);//Prob rem non lum aminoacid!
		copyState(&sn, s); //returns to state s variable sn
	}
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		if (info.dyeCanBeRemoved[i]) // A dye of that color can be removed.
		{
			pDA = pDyeAttached(sn, i);
			decN(remChar(&sn, i + '0'), i);
			setNewDyeSeqsIdxs(&sn, info.dyeSeqsIdxs[i] , info.dyeSeqsIdxsCount[i]); //Sets new dye sequences that correspond to state.
			if ((1-pDA)!=0)
				appendCandidateState(sn, (1-pDA) * info.pRemDye[i] * initProb);//Removed lum amino acid but wasnt attached!
			decK(&sn, i);
			if (pDA != 0)
				appendCandidateState(sn, pDA * info.pRemDye[i] * initProb);//Removed lum amino acid and was attached!
			copyState(&sn, s); //returns to state s variable sn
		}
	}

}

void Decoder::calcObservationProbs(float obs[N_COLORS])
{
	cw.getObsLogProbs(&auxObsLogProb, auxStates, obs);
	for (unsigned int i = 0; i < auxStates.size(); i++)
		auxStatesProb[i] = log(auxStatesProb[i]) + auxObsLogProb[i]; // Log prob of the state + log prob of observation.
	
}
void Decoder::keepBestStates(unsigned int t)
{
	vector<unsigned int> idxsSorted = argsortf(auxStatesProb);
	unsigned int currStateIdx;
	for (unsigned int i = 0; i < nBeam; i++)
	{
		currStateIdx = idxsSorted[i];
		copyState(&mostLikelyStates[t][i], auxStates[currStateIdx]);
		if (i == 0)
			mostLikelyStatesProbNorm[t][i] = auxStatesProb[currStateIdx] ;
		else
			mostLikelyStatesProbNorm[t][i] = auxStatesProb[currStateIdx] - mostLikelyStatesProbNorm[t][0];
	}
	mostLikelyStatesProbNorm[t][0] = 0; //Normalizes first value at the very end
}

pair<unsigned int, float> Decoder::getMostProbDyeSeqIdx()
{
	return cw.getMostProbDyeSeqIdx(mostLikelyStates[N_ED_CYC], mostLikelyStatesProbNorm[N_ED_CYC]);
}

void Decoder::appendCandidateState(State& s, float prob)
{
	bool found = false; //Pushes state if it was not here before, if it was adds the probability to transition to that state.
	if (prob < 0)
		cout << "This shouldnt happen";
	for (unsigned int i = 0; i < auxStates.size(); i++) 
	{
		if (isEqual(auxStates[i], s))
		{
			auxStatesProb[i] += prob; //If state was before, adds probabilities
			found = true;
			break;
		}
	}
	if (found == false)
	{
		auxStates.push_back(s);
		auxStatesProb.push_back(prob);
	}
}

void Decoder::clear()
{
	auxStates.clear();
	auxStatesProb.clear();
	cw.clear();
	/*for (unsigned int i = 0; i < N_FEATURES_PER_COL; i++) {
		mostLikelyStates[i].clear();
		mostLikelyStatesProbNorm[i].clear();
	}*/
}

//Useful functions outside of class
