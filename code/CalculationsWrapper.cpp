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
	dyeSeqsCounts = relCounts;
	reformatDyeSeqs(dyeSeqs, dyeSeqsIdx, relCounts);
	InfoForEdmanDegradation aux;
	for (unsigned int i = 0; i < nBeam; i++)
		infosEdman.push_back(aux);
	is.init(&dyeSeqsTogheter, &dyeSeqsStartIdxsInMem, &relProbs ,&dyeSeqsCounts, nBeam);
}

void CalculationsWrapper::reformatDyeSeqs(vector<string>& dyeSeqs, vector<unsigned int>& dyeSeqsIdx, vector<unsigned int>& relCounts)
{
	n_peptides = 0;
	unsigned long countChars = 0;
	//Reformatting data.
	
	for (unsigned int i = 0; i < relCounts.size(); i++)
	{
		relProbs.push_back(relCounts[i]);
		n_peptides += relCounts[i];
	}
		

	for (unsigned int i = 0; i < relProbs.size(); i++)
		relProbs[i] /= n_peptides;



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


float CalculationsWrapper::calcDyeLossProb(unsigned int Kinit[N_COLORS], unsigned int Kend[N_COLORS])
{
	float retVal = 0;
	unsigned int totalInitDyes = Kinit[0] + Kinit[1] + Kinit[2];
	float comb_factor = (float) (nChoosek(Kinit[0], Kend[0]) * nChoosek(Kinit[1], Kend[1]) * nChoosek(Kinit[2], Kend[2])); //Combinational factor of picking which dyes were not attached
	unsigned int dif = totalInitDyes - Kend[0] - Kend[1] - Kend[2];
	retVal=comb_factor * (float)pow(L, dif) * (float)pow((1 - L), totalInitDyes - dif);

	return retVal;
}

void CalculationsWrapper::getInfoForEdman(vector<State>& sV,unsigned int nStates) 
{
	string currDyeSeq;
	float currDyeProb;
	unsigned int dyeSeqIdx, dyeSeqLen;
	unsigned long dyeSeqStartInMemIdx;
	float totalProb = 0;
	bool probFound = false;
	char nextAA;
	for (unsigned int k = 0; k < nStates; k++) //For every state that we consider
	{
		State& s = sV[k];
		infosEdman[k].clear();//Clears variable
		for (unsigned int i = 0; i < s.dyeSeqsIdxsCount; i++) //Iterate over idxs
		{
			dyeSeqIdx = s.dyeSeqsIdxs[i]; //Index that represent a dye sequence
			dyeSeqStartInMemIdx = dyeSeqsStartIdxsInMem[dyeSeqIdx]; //Index that shows in our block memory (dyeSeqsTogheter) where our dye sequence starts 
			dyeSeqLen = dyeSeqsStartIdxsInMem[dyeSeqIdx + 1] - dyeSeqStartInMemIdx;
			currDyeProb = relProbs[dyeSeqIdx];
			if (dyeSeqLen - 1 > s.RCharCount) //In case there exists a possibility of removing an aminoacid
			{
				probFound = true;
				nextAA = dyeSeqsTogheter[dyeSeqStartInMemIdx + s.RCharCount];
				if (nextAA == '.')
				{
					infosEdman[k].nonLumCanBeRemoved = true;
					infosEdman[k].dyeSeqsIdxsDot[infosEdman[k].dyeSeqsIdxsCountDot++] = dyeSeqIdx;
					infosEdman[k].pRemNonLum += currDyeProb;
				}
				else
				{
					for (unsigned int j = 0; j < N_COLORS; j++)
					{
						if (nextAA == (j + '0'))
						{
							infosEdman[k].dyeCanBeRemoved[j] = true;
							infosEdman[k].dyeSeqsIdxs[j][infosEdman[k].dyeSeqsIdxsCount[j]++] = dyeSeqIdx;
							infosEdman[k].pRemDye[j] += currDyeProb;
						}
					}
				}
				totalProb += currDyeProb;
			}
		}
		if (probFound)
		{
			infosEdman[k].pRemNonLum /= totalProb;//If any of the probs then normalize probabilities
			for (unsigned int j = 0; j < N_COLORS; j++)
				infosEdman[k].pRemDye[j] /= totalProb;
		}
	}
	//if (abs((infoEdman.pRemNonLum + infoEdman.pRemDye[0] + infoEdman.pRemDye[1] + infoEdman.pRemDye[2])-1) > 0.1)
	//	cout << "This shouldnt happen";
	

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
				possibleOutputs[dyeSeqIdx] += finalStatesLogProbs[i] * dyeSeqsProbRelOut[d_idx] ;
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
		dyeSeqsProbRelOut[i] = 1/n_peptides;
		normFactor += relProbs[s.dyeSeqsIdxs[i]]; //
	}
	for (unsigned int i = 0; i < s.dyeSeqsIdxsCount; i++)
		dyeSeqsProbRelOut[i] /= normFactor;
}

void CalculationsWrapper::clear()
{
	is.clear();
	
}



