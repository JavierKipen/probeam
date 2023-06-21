#include "CalculationsWrapper.h"
#include <map>
#include <list>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include "StateFunctions.h"
#include <cmath>    

#define DYE_BY_MEM_IDX(i) (dyeSeqsTogheter[dyeSeqsStartIdxsInMem[(i)]])
#define L1NORM_MAX 5




unsigned int nChoosek(unsigned int n, unsigned int k);
void highestKValsInArray(float* array, unsigned long n, unsigned int k, unsigned long* outBestIdx);

CalculationsWrapper::CalculationsWrapper()
{
	l1normMaxDyeLoss = L1NORM_MAX; // Could be a parameter to run 
}


void CalculationsWrapper::init(vector<string>& dyeSeqs, vector<unsigned int>& dyeSeqsIdx, vector<unsigned int>& relCounts,unsigned int nBeam)
{
	l1normMaxDyeLoss = L1NORM_MAX; // Could be a parameter to run 
	//Mem reservation
	dyeSeqsTogheter.reserve(1000000);
	KDyeLoss.reserve(100);
	KProbsDyeLoss.reserve(100);
	relProbs.reserve(relCounts.size());
	dyeSeqsStartIdxsInMem.reserve(relCounts.size());
	reformatDyeSeqs(dyeSeqs, dyeSeqsIdx, relCounts);
	is.init(&dyeSeqsTogheter, &dyeSeqsStartIdxsInMem, &relProbs,nBeam);
}

void CalculationsWrapper::reformatDyeSeqs(vector<string>& dyeSeqs, vector<unsigned int>& dyeSeqsIdx, vector<unsigned int>& relCounts)
{
	float normFactorRel = 0;
	unsigned long countChars = 0;
	//Reformatting data.
	for (unsigned int i = 0; i < relCounts.size(); i++)
		normFactorRel += relCounts[i];
	for (unsigned int i = 0; i < relCounts.size(); i++)
		relProbs.push_back(((float)relCounts[i]) / normFactorRel);

	for (unsigned int i = 0; i < dyeSeqs.size(); i++)
	{
		dyeSeqsStartIdxsInMem.push_back(countChars);
		copy(dyeSeqs[i].begin(), dyeSeqs[i].end(), back_inserter(dyeSeqsTogheter));
		dyeSeqsTogheter.push_back('\n');
		countChars += dyeSeqs[i].size() + 1;
	}
	dyeSeqsStartIdxsInMem.push_back(countChars); //Final ID (no state should use it, but simplifies some calculations)
	dyeSeqsIdxOUT = dyeSeqsIdx;
}

void CalculationsWrapper::getFirstMostLikelyStates(vector<State>* outMostLikely, vector<float>* outProbsNorm,  float obs[N_COLORS])
{
	is.getFirstMostLikelyStates(outMostLikely, outProbsNorm, obs);
}

string CalculationsWrapper::internalDyeSeqIdToStr(unsigned int id)
{
	unsigned long startInMem = dyeSeqsStartIdxsInMem[id];
	string retVal;
	while (startInMem < dyeSeqsStartIdxsInMem[id + 1]-1)
		retVal.push_back(dyeSeqsTogheter[startInMem++]);
	return retVal;
}


void CalculationsWrapper::obtainKDyeLoss(State& s, float obs[N_COLORS])
{
	unsigned int count = 0;
	float comb_factor;
	unsigned int dif;
	unsigned int totalInitDyes = s.K[0] + s.K[1] + s.K[2];

	float L1dif_i, L1dif_j, L1dif_k;

	KDyeLoss.clear();
	KProbsDyeLoss.clear();

	for (unsigned int i = 0; i < s.K[0] + 1; i++) { //Loops with all possible Ks that could be
		L1dif_i = abs((obs[0] / MU) - i);
		for (unsigned int j = 0; j < s.K[1] + 1; j++) {
			L1dif_j = abs((obs[1] / MU) - j);
			for (unsigned int k = 0; k < s.K[2] + 1; k++) {
				L1dif_k = abs((obs[2] / MU) - k);
				if ((L1dif_i + L1dif_j + L1dif_k) < l1normMaxDyeLoss) //This K can give an observation probability which is not negligible
				{
					comb_factor = (float)nChoosek(s.K[0], i) * nChoosek(s.K[1], j) * nChoosek(s.K[2], k); //Combinational factor of picking which dyes were not attached
					dif = totalInitDyes - i - j - k;
					array<unsigned int, N_COLORS> aux = { i, j, k };
					KDyeLoss.push_back(aux);
					KProbsDyeLoss.push_back((float)comb_factor * (float)pow(L, dif) * (float)pow((1 - L), totalInitDyes - dif));//Probability of initial state
				}
			}
		}
	}
}

void CalculationsWrapper::getInfoForEdman(State& s)
{
	string currDyeSeq;
	float currDyeProb;
	unsigned int dyeSeqIdx,dyeSeqLen;
	unsigned long dyeSeqStartInMemIdx;
	float totalProb=0;
	bool probFound = false;
	char nextAA;
	infoEdman.clear();//Clears variable
	for (unsigned int i = 0; i < s.dyeSeqsIdxsCount; i++) //Iterate over idxs
	{
		dyeSeqIdx = s.dyeSeqsIdxs[i]; //Index that represent a dye sequence
		dyeSeqStartInMemIdx = dyeSeqsStartIdxsInMem[dyeSeqIdx]; //Index that shows in our block memory (dyeSeqsTogheter) where our dye sequence starts 
		dyeSeqLen = dyeSeqsStartIdxsInMem[dyeSeqIdx + 1] - dyeSeqStartInMemIdx;
		currDyeProb = relProbs[dyeSeqIdx];
		if (dyeSeqLen-1 > s.RCharCount) //In case there exists a possibility of removing an aminoacid
		{
			probFound = true;
			nextAA = dyeSeqsTogheter[dyeSeqStartInMemIdx + s.RCharCount];
			if (nextAA == '.')
			{
				infoEdman.nonLumCanBeRemoved = true;
				infoEdman.dyeSeqsIdxsDot[infoEdman.dyeSeqsIdxsCountDot++] = dyeSeqIdx;
				infoEdman.pRemNonLum += currDyeProb;
			}
			else
			{
				for (unsigned int j = 0; j < N_COLORS; j++)
				{
					if (nextAA == (j + '0'))
					{
						infoEdman.dyeCanBeRemoved[j] = true;
						infoEdman.dyeSeqsIdxs[j][infoEdman.dyeSeqsIdxsCount[j]++]= dyeSeqIdx;
						infoEdman.pRemDye[j] += currDyeProb;
					}
				}
			}
			totalProb += currDyeProb;
		}
	}
	if (probFound)
	{
		infoEdman.pRemNonLum /= totalProb;//If any of the probs then normalize probabilities
		for (unsigned int j = 0; j < N_COLORS; j++)
			infoEdman.pRemDye[j] /= totalProb;
	}
	//if (abs((infoEdman.pRemNonLum + infoEdman.pRemDye[0] + infoEdman.pRemDye[1] + infoEdman.pRemDye[2])-1) > 0.1)
	//	cout << "This shouldnt happen";
	

}

void CalculationsWrapper::getObsLogProbs(vector<float>* outLogProbs, vector<State>& auxStates, float obs[N_COLORS])
{
	float stds_loc[N_COLORS];
	float means_loc[N_COLORS];
	float auxLogProbObs;

	for (unsigned int i = 0; i < auxStates.size(); i++)
	{
		State& Si = auxStates[i];
		auxLogProbObs = 0;
		means_loc[0] = getMeanState(Si, 0); means_loc[1] = getMeanState(Si, 1); means_loc[2] = getMeanState(Si, 2);
		stds_loc[0] = getStdState(Si, 0); stds_loc[1] = getStdState(Si, 1); stds_loc[2] = getStdState(Si, 2);
		for (unsigned long j = 0; j < N_COLORS; j++)
			auxLogProbObs -= (float)pow((obs[j] - means_loc[j]) / stds_loc[j], 2);
		auxLogProbObs /= 2;
		auxLogProbObs -= (float)log(LOGTERM_CONST * stds_loc[0] * stds_loc[1] * stds_loc[2]);
		outLogProbs->push_back(auxLogProbObs); //Pushes back logprob of obs
	}
}

pair<unsigned int, float> CalculationsWrapper::getMostProbDyeSeqIdx(vector<State>& finalStates, vector<float>& finalStatesLogProbs)
{
	map<unsigned int, float> possibleOutputs;
	pair<unsigned int, float> output;
	float normFactor = 0;
	unsigned int nBeam = finalStatesLogProbs.size();
	unsigned int mostLikelyOut;
	float mostLikelyOutP = -1;
	//Normalizing last state probabilities
	for (unsigned int i = 0; i < nBeam; i++)
	{
		finalStatesLogProbs[i] = exp(finalStatesLogProbs[i]);
		normFactor += finalStatesLogProbs[i];
	}
	for (unsigned int i = 0; i < nBeam; i++)
		finalStatesLogProbs[i] /= normFactor; //Normalizes probabilities of states
	//State probabilities to dye seq probabilities
	for (unsigned int i = 0; i < nBeam; i++)
	{
		State& s = finalStates[i]; //Ith most probable state at the end of the observation
		getRelProbs(s); // Computes the rel prob of each dye sequence and saves it in dyeSeqsProbRelOut;
		for (unsigned int d_idx = 0; d_idx < s.dyeSeqsIdxsCount; d_idx++)
		{
			unsigned int dyeSeqIdx = s.dyeSeqsIdxs[d_idx];
			auto it = possibleOutputs.find(dyeSeqIdx);
			if (it == possibleOutputs.end()) //Idx of dye seq not found
				possibleOutputs[dyeSeqIdx] = finalStatesLogProbs[i] * dyeSeqsProbRelOut[d_idx];
			else
				possibleOutputs[dyeSeqIdx] += finalStatesLogProbs[i] * dyeSeqsProbRelOut[d_idx];
		}

	}
	//Picking most likely
	for (auto it = possibleOutputs.begin(); it != possibleOutputs.end(); it++)
	{
		if (it->second > mostLikelyOutP)
		{
			mostLikelyOut = it->first;
			mostLikelyOutP = it->second;
		}
	}
	string picked = internalDyeSeqIdToStr(mostLikelyOut);
	output.first = dyeSeqsIdxOUT[mostLikelyOut];
	output.second = mostLikelyOutP;
	return output;
}

void  CalculationsWrapper::getRelProbs(State& s)
{
	float normFactor=0;
	for (unsigned int i = 0; i < s.dyeSeqsIdxsCount; i++)
	{
		dyeSeqsProbRelOut[i] = relProbs[s.dyeSeqsIdxs[i]];
		normFactor += dyeSeqsProbRelOut[i];
	}
	for (unsigned int i = 0; i < s.dyeSeqsIdxsCount; i++)
		dyeSeqsProbRelOut[i] /= normFactor;
}

void CalculationsWrapper::clear()
{
	is.clear();
}



