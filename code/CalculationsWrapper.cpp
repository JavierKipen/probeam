#include "CalculationsWrapper.h"
#include <map>
#include <list>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include "StateFunctions.h"
#include <cmath>    
#include <chrono>

#define DYE_BY_MEM_IDX(i) (dyeSeqsTogheter[dyeSeqsStartIdxsInMem[(i)]])
#define L1NORM_MAX 5




unsigned int nChoosek(unsigned int n, unsigned int k);
void highestKValsInArray(float* array, unsigned long n, unsigned int k, unsigned long* outBestIdx);
//vector<unsigned int> argsortf(const vector<float>& vf);
bool customSort(const pair<unsigned int, string> firstElem, const pair<unsigned int, string> secondElem);

CalculationsWrapper::CalculationsWrapper()
{
	l1normMaxDyeLoss = L1NORM_MAX; // Could be a parameter to run 
}


void CalculationsWrapper::init(vector<string>& dyeSeqs, vector<unsigned int>& dyeSeqsIdx, vector<unsigned int>& relCounts,unsigned int nBeam)
{
	l1normMaxDyeLoss = L1NORM_MAX; // Could be a parameter to run 
	vector<float> vect(dyeSeqsIdx.size(), 0);
	dyeSeqsOutProb = vect;
	//Mem reservation
	dyeSeqsTogheter.reserve(1000000);
	KDyeLoss.reserve(100);
	KProbsDyeLoss.reserve(100);
	dyeSeqsOut.reserve(20000); //Vector to calculate the final peptide probs
	relProbs.reserve(relCounts.size());
	dyeSeqsStartIdxsInMem.reserve(relCounts.size());
	dyeSeqsCounts = relCounts;

	reOrderDyeSeqs(dyeSeqs, dyeSeqsIdx, relCounts);
	//reformatDyeSeqs(dyeSeqs, dyeSeqsIdx, relCounts);
	InfoForEdmanDegradation aux;
	for (unsigned int i = 0; i < nBeam; i++)
		infosEdman.push_back(aux);
	is.init(&dyeSeqsTogheter, &dyeSeqsStartIdxsInMem, &relProbs ,&dyeSeqsCounts, nBeam);

#ifdef ESTIMATE_IFED_TIMES
	initMapTimeIFED();
#endif // ESTIMATE_IFED_TIMES

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
		#ifdef ESTIMATE_IFED_TIMES
			auto start = chrono::high_resolution_clock::now();
		#endif // ESTIMATE_IFED_TIMES
		totalProb = 0;
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
	#ifdef ESTIMATE_IFED_TIMES
		auto stop = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
		IFEDTimeOfState[getStateOriginalN(s)] += ((float)duration.count()) * 1e-9;
	#endif // ESTIMATE_IFED_TIMES
	}
	//if (abs((infoEdman.pRemNonLum + infoEdman.pRemDye[0] + infoEdman.pRemDye[1] + infoEdman.pRemDye[2])-1) > 0.1)
	//	cout << "This shouldnt happen";
	

}


pair<unsigned int, float> CalculationsWrapper::getMostProbDyeSeqIdx(vector<State>& finalStates, vector<float>& finalStatesLogProbs,unsigned int currNStates)
{
	pair<unsigned int, float> output;
	float normFactor = 0;
	unsigned int nBeam = currNStates;
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
			unsigned int outDyesIdx = 0;
			dyeSeqsOut.push_back(dyeSeqIdx);
			dyeSeqsOutProb[dyeSeqIdx] += (finalStatesLogProbs[i] * dyeSeqsProbRelOut[d_idx]);
		}

	}
	//Remove repeated sequences
	vector<unsigned int>::iterator ip;
	ip = unique(dyeSeqsOut.begin(), dyeSeqsOut.end());
	dyeSeqsOut.resize(distance(dyeSeqsOut.begin(), ip));

	//Picking most likely
	for (unsigned int i=0;i< dyeSeqsOut.size();i++)
	{
		if (dyeSeqsOutProb[dyeSeqsOut[i]] > mostLikelyOutP)
		{
			mostLikelyOut = dyeSeqsOut[i];
			mostLikelyOutP = dyeSeqsOutProb[dyeSeqsOut[i]];
		}
		dyeSeqsOutProb[dyeSeqsOut[i]] = 0; //Resets output var
	}
	//string picked = internalDyeSeqIdToStr(mostLikelyOut);
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

void  CalculationsWrapper::initMapTimeIFED()
{
	array<unsigned int, N_COLORS> auxN;
	for (auto it : is.initIdealStates)
	{
		for (unsigned int i = 0; i < N_COLORS; i++)
			auxN[i] = it.N[i];
		IFEDTimeOfState[auxN] = 0;
	}
}

array<unsigned int, N_COLORS> CalculationsWrapper::getStateOriginalN(State& s)
{
	array<unsigned int, N_COLORS> origN;
	for (unsigned int i = 0; i < N_COLORS; i++)
		origN[i] = s.N[i];
	for (unsigned int i = 0; i < s.RCharCount; i++)
	{
		if (s.R[i] != '.')
		{
			origN[(s.R[i] - '0')]++;
		}
	}
	return origN;
}

void CalculationsWrapper::orderCompTimesIFED()
{
	vector<array<unsigned int, N_COLORS>> IdealStateNs;
	vector<float> PercTimeSpent;
	float normVal = 0;
	for (auto it = IFEDTimeOfState.begin(); it != IFEDTimeOfState.end(); ++it)
	{
		IdealStateNs.push_back(it->first);
		PercTimeSpent.push_back(it->second);
		normVal += it->second;
	}
	for (unsigned int i = 0; i < PercTimeSpent.size(); i++)
		PercTimeSpent[i] /= normVal;
	vector<unsigned int> IdxSort = argsortf(PercTimeSpent);

	vector<float> PercTimeSpentOrdered= PercTimeSpent;
	vector<array<unsigned int, N_COLORS>> IdealStateNsOrdered = IdealStateNs;
	for (unsigned int i = 0; i < PercTimeSpent.size(); i++)
	{
		PercTimeSpentOrdered[i] = PercTimeSpent[IdxSort[i]];
		IdealStateNsOrdered[i] = IdealStateNs[IdxSort[i]];
	}

}


void CalculationsWrapper::clear()
{
	is.clear();
	dyeSeqsOut.clear(); //Vector to calculate the final peptide probs

}


//Reordering of dye sequences in memory 

void CalculationsWrapper::reOrderDyeSeqs(vector<string>& dyeSeqs, vector<unsigned int>& peptideIdx, vector<unsigned int>& relCounts)
{
	map<array<unsigned int, N_COLORS>, list<pair<unsigned int,string>>> groupingDyeSeqs;
	n_peptides = 0;
	unsigned long countChars = 0;

	for (unsigned int i = 0; i < peptideIdx.size(); i++)
	{
		pair<unsigned int, string> aux(i, dyeSeqs[i]);
		groupingDyeSeqs[getNDyeSeq(dyeSeqs[i])].push_back(aux);
	}
		
	//Here N could be ordered in a specific way too.
	for (auto Nit : groupingDyeSeqs)
	{
		list<pair<unsigned int, string>> dyeSeqsForGivenN = Nit.second;
		dyeSeqsForGivenN.sort(customSort);
		for (auto it = dyeSeqsForGivenN.begin(); it != dyeSeqsForGivenN.end(); ++it)
		{
			string s = it->second;
			unsigned int origDyeSeqIdx = it->first;
			relProbs.push_back(relCounts[origDyeSeqIdx]);
			n_peptides += relCounts[origDyeSeqIdx];
			dyeSeqsStartIdxsInMem.push_back(countChars); //Start of the new index
			copy(s.begin(), s.end(), back_inserter(dyeSeqsTogheter));
			dyeSeqsTogheter.push_back('\n');
			countChars += s.size() + 1;
			dyeSeqsIdxOUT.push_back(peptideIdx[origDyeSeqIdx]);
		}
	}
	dyeSeqsStartIdxsInMem.push_back(countChars); //Final ID (no state should use it, but simplifies some calculations)
	for (unsigned int i = 0; i < relProbs.size(); i++)
		relProbs[i] /= (float)n_peptides;
}

array<unsigned int, N_COLORS> CalculationsWrapper::getNDyeSeq(string &dyeSeq)
{
	array<unsigned int, N_COLORS> retVal;
	for(unsigned int i=0;i<N_COLORS;i++)
		retVal[i]=count(dyeSeq.begin(), dyeSeq.end(), '0'+i);
	return retVal;
}

bool customSort(const pair<unsigned int, string> firstElem, const pair<unsigned int, string> secondElem)
{
	return firstElem.second< secondElem.second;
}

/*

vector<unsigned int> argsortf(const vector<float>& vf) { //Argsort for a float vector. based on https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
	vector<unsigned int> idx(vf.size());// initialize original index locations
	iota(idx.begin(), idx.end(), 0);
	stable_sort(idx.begin(), idx.end(),
		[&vf](unsigned int i1, unsigned int i2) {return vf[i1] > vf[i2]; });
	return idx;
}*/