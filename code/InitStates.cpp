#include "InitStates.h"
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <cmath>    
#include <array>

//#define WHOLE_CALC_INITIAL //When this is defined, all observations are calculated for the initial states
#define BIG_FLOAT 1E20 // Big float to intiialize values ( FLT_MAX did not compile with g++ idk why)

void getPosNextKValues(vector<array<unsigned int, N_COLORS>>& possibleNextKs, unsigned int sKOrig[N_COLORS], unsigned int sKBest[N_COLORS], unsigned int N[N_COLORS]);



InitStates::InitStates()
{
	bestStatesInfo.reserve(100); //Nbeam usually lower than this
	logL = log(L);
	logOneMinusL = log(1-L);
}


void InitStates::init(vector<char> *dyeSeqsTogheterO, vector<unsigned long>* dyeSeqsStartIdxsInMemO, vector<float>* relProbsO, vector<unsigned int>* dyeSeqsCounts, unsigned int nBeamIn)
{
	nDyeSeqs = dyeSeqsStartIdxsInMemO->size();
	dyeSeqsTogheter = dyeSeqsTogheterO;
	dyeSeqsStartIdxsInMem = dyeSeqsStartIdxsInMemO;
	relProbs = relProbsO;
	nBeam = nBeamIn;
	dyeSeqsCountsp = dyeSeqsCounts;

	createInitIdealStates();
#ifdef WHOLE_CALC_INITIAL
	extendToInitStates();
	precomputeForHardcoreInitStates();
#else
	getMaxNi();
	generatePosKsInit();
	possibleNextKs.reserve(200);
	possibleNextKsAux.reserve(200);
#endif 
}

void InitStates::createInitIdealStates()
{
	unsigned int currN[N_COLORS];
	unsigned int dyeSeqCount = 0;
	long charCount = 0;
	vector<char>& auxDyeSeqTog = *dyeSeqsTogheter;
	while (dyeSeqCount < nDyeSeqs-1)
	{
		for (unsigned int c = 0; c < N_COLORS; c++)
			currN[c] = 0; //resets color counter
		while (auxDyeSeqTog[charCount] != '\n')
		{
			for (unsigned int c = 0; c < N_COLORS; c++)
			{
				if (auxDyeSeqTog[charCount] == c + '0')
					currN[c]++;
			}
			charCount++;
		}
		appendForIdInSt(currN, dyeSeqCount, (*relProbs)[dyeSeqCount] );
		dyeSeqCount++; charCount++;
	}
	for (unsigned int i = 0; i < initIdealStatesProb.size(); i++)
		initIdealStatesLogProb.push_back(log(initIdealStatesProb[i]));
}

void InitStates::appendForIdInSt(unsigned int N[N_COLORS], unsigned int dyeSeqIdx, float prob)
{
	bool found = false; //Pushes state if it was not here before, if it was adds the probability to transition to that state.
	State auxState;
	auxState.dyeSeqsIdxsCount = 1; //Whenever pushing a state, only 1 dye seq is present.
	auxState.RCharCount = 0; //No removal on init ideal states
	for (unsigned int i = 0; i < initIdealStates.size(); i++)
	{
		if (isNEqual(initIdealStates[i], N))
		{
			initIdealStatesProb[i] += prob; //If state was before, adds probabilities
			initIdealStates[i].dyeSeqsIdxs[initIdealStates[i].dyeSeqsIdxsCount++] = dyeSeqIdx; //Adds dye seq
			found = true;
			break;
		}
	}
	if (found == false)
	{
		setK(setN(&auxState, N), N);
		auxState.dyeSeqsIdxs[0] = dyeSeqIdx;
		initIdealStatesProb.push_back(prob);
		initIdealStates.push_back(auxState);
	}
}

void InitStates::getFirstMostLikelyStates(vector<State>* outMostLikely, vector<float>* outProbsNorm,  float obs[N_COLORS])
{
#ifdef WHOLE_CALC_INITIAL
	hardcoreMostLikelyInitStates(outMostLikely, outProbsNorm, nBeam, obs);
#else
	smartMostLikelyInitStates(outMostLikely, outProbsNorm, obs);
#endif 
	//else the others.
}



void InitStates::clear()
{
	#ifdef WHOLE_CALC_INITIAL
		initStatesLogProbHCopy = initStatesLogProbH;
	#else
		bestStatesInfo.clear();
		obsProbMap.clear();
	#endif 
	
}

//For smart init states (minimizing computing time and memory usage)
void InitStates::getMaxNi()
{

	for (unsigned int i = 0; i < N_COLORS; i++)
		maxNis[i] = 0;
	for (unsigned int i = 0; i < initIdealStates.size(); i++)
		for (unsigned int j = 0; j < N_COLORS; j++)
			if (initIdealStates[i].N[j] > maxNis[j])
				maxNis[j] = initIdealStates[i].N[j];
	for (unsigned int j = 0; j < N_COLORS; j++)
		shiftsMin[j] = (unsigned int)floor(log2(maxNis[j])) + 1;

}

unsigned int InitStates::NorKToUnsInt(unsigned int NorK[N_COLORS])
{
	unsigned int retVal = NorK[2];
	retVal = (retVal << shiftsMin[1]) + NorK[1];
	retVal = (retVal << shiftsMin[0]) + NorK[0];
	return retVal;
}
void InitStates::KuintToVector(unsigned int NorK[N_COLORS], unsigned int baseK)
{
	NorK[0] = (baseK % (1 <<shiftsMin[0]));
	NorK[1] = (baseK>> shiftsMin[0]) % (1 << shiftsMin[1]);
	NorK[2] = baseK >> (shiftsMin[0] + shiftsMin[1]);
}

void InitStates::generatePosKsInit()
{
	float comb_factor, aux_prob, dif, totalInitDyes;
	map<unsigned int, float> auxMapProbKs;
	unsigned int auxK[N_COLORS];
	PdyeMissInit.reserve(initIdealStates.size());
	for (auto& s : initIdealStates)
	{
		auxMapProbKs.clear();
		for (unsigned int i = 0; i < s.N[0] + 1; i++) { //Loops with all possible Ks that could be
			for (unsigned int j = 0; j < s.N[1] + 1; j++) {
				for (unsigned int k = 0; k < s.N[2] + 1; k++) {
					if (i != 0 || j != 0 || k != 0) //Does not consider 000 states. 
					{
						auxK[0] = i; auxK[1] = j; auxK[2] = k;
						comb_factor = (float)nChoosek(s.N[0], i) * (float)nChoosek(s.N[1], j) * (float)nChoosek(s.N[2], k); //Combinational factor of picking which dyes were not attached
						totalInitDyes = s.N[0] + s.N[1] + s.N[2];
						dif = totalInitDyes - i - j - k;
						aux_prob = comb_factor * pow(M, dif) * pow((1 - M), totalInitDyes - dif); //Relative property within the possible Ks
						auxMapProbKs[NorKToUnsInt(auxK)]=(float) log(aux_prob); //Probability of initial state
					}
				}
			}
		}
		PdyeMissInit.push_back(auxMapProbKs);
		//To obtain the max prob distance
		float bestDiffs = 0;
		for (unsigned int c = 0; c < N_COLORS; c++) { //For each color I search the maximum increment in probability that could happen
			for (unsigned int i = 1; i < s.N[c] + 1; i++) {
				float obsKPrev = nChoosek(s.N[c], i - 1) * pow(M, s.N[c] - (i - 1)) * pow((1 - M), (i-1));
				float obsKNext = nChoosek(s.N[c],   i  ) * pow(M,     s.N[c] - i  ) * pow((1 - M),   i  );
				if ((obsKNext - obsKPrev) > bestDiffs)
					bestDiffs = (obsKNext - obsKPrev);
			}
		}
		maxIncDyeMiss.push_back(bestDiffs);
	}
}

vector<unsigned int> InitStates::orderByClosestStates(float obs[N_COLORS])
{
	vector<float> distSq(initStates.size());
	vector<unsigned int> retVal;
	float obs_norm[N_COLORS];
	float aux_dist;
	for (unsigned int i=0; i < N_COLORS; i++)
		obs_norm[i] = obs[i] / MU;
	for (unsigned int i=0; i < initIdealStates.size(); i++)
	{
		aux_dist = 0;
		for (unsigned int j=0; j < N_COLORS; j++)
			aux_dist += SQUARE((obs_norm[j] - initIdealStates[i].N[j]));
		distSq.push_back(aux_dist);
	}
	retVal = argsortf(distSq);
	reverse(retVal.begin(), retVal.end());
	return retVal;
}

void InitStates::smartMostLikelyInitStates(vector<State>* outMostLikely, vector<float>* outProbsNorm,  float obs[N_COLORS])
{
	bestStatesIdx = orderByClosestStates(obs); //Orders states by closeness to the obs.
	getBestK(obs); //Gets more likely Ks of observations.
	preLoadBestStates(obs); //bestStatesInfo now will contain the guess of 10 most likely states
	smartSearchStates(obs); //Looks wisely through all states and K to see which are the most likely states
	formatOutput(outMostLikely, outProbsNorm); // From the best states in this class specific format obtains the states for the decoder format
}

void InitStates::formatOutput(vector<State>* outMostLikely, vector<float>* outProbsNorm)
{
	State aux;
	unsigned int count = 0;
	float normFloat = bestStatesInfo[0].logProb;
	for (auto& isi : bestStatesInfo) //Init States Info
	{
		copyState(&aux, initIdealStates[isi.stateIndex]);
		KuintToVector(aux.K, isi.K);
		(*outMostLikely)[count]= aux;
		(*outProbsNorm)[count++]=(isi.logProb - normFloat);
	}
}
void InitStates::getBestK(float obs[N_COLORS])
{
	unsigned int auxKi;
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		auxKi = (int(obs[i]/MU) < 0) ? 0 : int(obs[i]/MU); //Lowest possible value is floor or 0.
		KBest[i]=(probObsSingleColor(obs[i], auxKi) > probObsSingleColor(obs[i], auxKi + 1)) ? auxKi : auxKi+1; //Keeps the most likely Ki
	}

}


void InitStates::preLoadBestStates(float obs[N_COLORS])
{
	vector<float> logProbs;
	vector<InitStateInfo> bestStatesInfoAux; //To keep track of the states but unordered.
	float logProb;
	InitStateInfo aux;
	unsigned int sKBest[N_COLORS];
	for (unsigned int i = 0; i < nBeam; i++)
	{
		unsigned int sIdx = bestStatesIdx[i]; //Index of the state to work with
		getBestKState(sKBest, initIdealStates[sIdx].N); // We get the best K for the state (Kbest if N is higher than it for every value).
		aux.K = NorKToUnsInt(sKBest); //Converts the format of array to single uint
		aux.stateIndex = sIdx;
		logProb = initIdealStatesLogProb[sIdx] + PdyeMissInit[sIdx][aux.K] + getObsLogProb(obs, aux.K); //log( S Init prob * PDyeMiss * Pobs );
		aux.logProb = logProb;
		bestStatesInfoAux.push_back(aux);
		logProbs.push_back(logProb);
	}
	vector<unsigned int> order = argsortf(logProbs); // with the order now we push to best states
	for (unsigned int i = 0; i < order.size(); i++) 
		bestStatesInfo.push_back(bestStatesInfoAux[order[i]]); //Pushes the ordered first guess of states

}

void InitStates::smartSearchStates( float obs[N_COLORS])
{
	InitStateInfo aux;
	unsigned int sKBest[N_COLORS];
	for (unsigned int i = 0; i < bestStatesIdx.size(); i++)
	{
		unsigned int sIdx = bestStatesIdx[i];
		getBestKState(sKBest, initIdealStates[sIdx].N); // We get the best K for the state (Kbest if N is higher than it for every value).
		aux.stateIndex = sIdx;
		aux.K = NorKToUnsInt(sKBest); //Converts the format of array to single uint
		if (i > nBeam) //For the first nBeam its not needed to check the best Kstate, since it was used to fill the buffer
		{
			float stateNObsLog= initIdealStatesLogProb[sIdx] + getObsLogProb(obs, aux.K); //P(S) * Pobs with Kopt
			if (stateNObsLog > bestStatesInfo[nBeam - 1].logProb) //We only check for neighbours if this prob is higher than the worst state of the list
			{
				//aux.logProb = stateNObsLog + PdyeMissInit[sIdx][aux.K]; //LogProb of the state
				aux.logProb = stateNObsLog + calcDyeLossProbLog(initIdealStates[sIdx].N, sKBest); //LogProb of the state
				if (aux.logProb > bestStatesInfo[nBeam - 1].logProb) //If the probability is higher than the last one
					pushInitStateInfo(aux); //Pushes the obtained value
				exploreOtherKs(aux, sKBest, obs); //Explores other Ks
			}
		}
		else
			exploreOtherKs(aux, sKBest, obs); //Explores other Ks
		
	}
}
float InitStates::getObsLogProb(float obs[N_COLORS], unsigned int K_uint)
{
	float retVal;
	unsigned int Kaux[N_COLORS];
	/*
	KuintToVector(Kaux, K_uint);
	retVal = logpObs(Kaux, obs);*/

	if (obsProbMap.find(K_uint) == obsProbMap.end()) {
		
		KuintToVector(Kaux, K_uint);
		retVal = logpObs(Kaux, obs);
		obsProbMap[K_uint] = retVal;
	}
	else {
		retVal = obsProbMap[K_uint];
	}
	return retVal;
}
void InitStates::getBestKState(unsigned int K[N_COLORS], unsigned int N[N_COLORS])
{
	for (unsigned int i = 0; i < N_COLORS; i++)
		K[i] = (N[i] > KBest[i]) ? KBest[i] : N[i]; //Keeps the most likely Ki
}

void InitStates::pushInitStateInfo(InitStateInfo& aux)
{
	unsigned int indexPos;
	for (indexPos = nBeam - 1; indexPos > 0; indexPos--)
		if (aux.logProb < bestStatesInfo[indexPos].logProb)
			break;
	indexPos++;
	if (indexPos == 1 && aux.logProb > bestStatesInfo[0].logProb)
		indexPos = 0;
	for (unsigned int i = nBeam - 1; i > indexPos; i--)
		bestStatesInfo[i] = bestStatesInfo[i - 1];
	bestStatesInfo[indexPos] = aux;
}

void InitStates::exploreOtherKs(InitStateInfo& stateInfo,unsigned int sKBest[N_COLORS], float obs[N_COLORS]) //breadth first search
{
	possibleNextKs.clear(); //For exploration of other Ks!
	possibleNextKsAux.clear();
	
	bool finishedSearch = false;
	getPosNextKValues(possibleNextKs, sKBest, sKBest, initIdealStates[stateInfo.stateIndex].N);
	InitStateInfo aux = stateInfo;
	while (!possibleNextKs.empty())
	{
		for (auto& newK : possibleNextKs)
		{
			aux.K = NorKToUnsInt(newK.data());
			float stateNObsLog = initIdealStatesLogProb[aux.stateIndex] + getObsLogProb(obs, aux.K); //P(S) * Pobs with Kopt
			if (stateNObsLog > bestStatesInfo[nBeam - 1].logProb) //If it could be that with this K a neighbour can have better prob than what we got
			{
				getPosNextKValues(possibleNextKsAux, newK.data(), sKBest, initIdealStates[aux.stateIndex].N);
				//aux.logProb = stateNObsLog + PdyeMissInit[aux.stateIndex][aux.K]; //LogProb of the state
				aux.logProb = stateNObsLog + calcDyeLossProbLog(initIdealStates[aux.stateIndex].N,newK.data()); //LogProb of the state
				if (aux.logProb > bestStatesInfo[nBeam - 1].logProb) //If the probability is higher than the last one
					pushInitStateInfo(aux); //Pushes the obtained value
				
			}
		}
		possibleNextKs.clear();
		vector<array<unsigned int, N_COLORS>>::iterator it = unique(possibleNextKsAux.begin(), possibleNextKsAux.end()); //Items can be repeated!
		for (auto auxIt = possibleNextKsAux.begin(); auxIt != it; ++auxIt)
			possibleNextKs.push_back(*auxIt);
		possibleNextKsAux.clear();
		//vector<array<unsigned int, N_COLORS>>::iterator it = unique(possibleNextKsAux.begin(), possibleNextKsAux.end()); //Items can be repeated!
		//possibleNextKsAux.resize(distance(possibleNextKsAux.begin(), it)); // 10 20 30 20 10
		//possibleNextKs = possibleNextKsAux;
		//possibleNextKsAux.clear();
	}
}

void getPosNextKValues(vector<array<unsigned int, N_COLORS>> &possibleNextKs, unsigned int sKOrig[N_COLORS], unsigned int sKBest[N_COLORS], unsigned int N[N_COLORS])
{
	array<unsigned int, N_COLORS> auxK;
	unsigned int initDist = KDist(sKOrig, sKBest);
	for (unsigned int i = 0; i < N_COLORS; i++)
		auxK[i] = sKOrig[i]; //Copy the origin K
	for (unsigned int i = 0; i < N_COLORS; i++) //For every color we see if when we get 1 fluorophore distance should we consider or not.
	{
		
		if (sKOrig[i] < N[i]) //Whenever by adding 1 we dont reach the limit of N
		{
			auxK[i] = auxK[i]+ 1; //If it was picked, we will save the Ki with a +1 value
			if(!isNEqual(auxK.data(), sKBest) && KDist(auxK.data(), sKBest)>initDist) //Push only of its not the base one and it is further away than before.
				possibleNextKs.push_back(auxK);
			auxK[i] = auxK[i] - 1;//To restore to same value
		}
		
		if (sKOrig[i] > 0) //We can decrease a K only if it is higher than 0.
		{
			auxK[i] = auxK[i] - 1;
			if (!isNEqual(auxK.data(), sKBest) && KDist(auxK.data(), sKBest) > initDist) //Push only of its not the base one and it is further away than before.
				possibleNextKs.push_back(auxK);
			auxK[i] = auxK[i] + 1;//To restore to same value
		}
	}
}

//For Hardcore init states (checking every single probability with every posible state.

void InitStates::extendToInitStates()
{
	unsigned long auxCount, initIdealStateCount;
	unsigned int dif, totalInitDyes;
	double aux_prob, comb_factor;
	unsigned int auxKValue[N_COLORS];
	auxCount = 0; initIdealStateCount = 0;
	State auxS;
	for (auto& is : initIdealStates)
	{
		copyState(&auxS, is);
		for (unsigned int i = 0; i < is.N[0] + 1; i++) { //Loops with all possible Ks that could be
			for (unsigned int j = 0; j < is.N[1] + 1; j++) {
				for (unsigned int k = 0; k < is.N[2] + 1; k++) {
					if (i != 0 || j != 0 || k != 0) //Does not consider 000 states. 
					{
						auxKValue[0] = i; auxKValue[1] = j; auxKValue[2] = k;
						initStates.push_back(*setK(&auxS, auxKValue));
						comb_factor = (double)nChoosek(is.N[0], i) * (double)nChoosek(is.N[1], j) * (double)nChoosek(is.N[2], k); //Combinational factor of picking which dyes were not attached
						totalInitDyes = is.N[0] + is.N[1] + is.N[2];
						dif = totalInitDyes - i - j - k;
						aux_prob = comb_factor * pow(M, dif) * pow((1 - M), totalInitDyes - dif); //Relative property within the possible Ks
						initStatesLogProbH.push_back((float)log(aux_prob * initIdealStatesProb[initIdealStateCount])); //Probability of initial state
					}
				}
			}
		}
		initIdealStateCount++;
	}
	initStatesLogProbHCopy = initStatesLogProbH;
}

void InitStates::hardcoreMostLikelyInitStates(vector<State>* outMostLikely, vector<float>* outProbsNorm, unsigned int nBeam, float obs[N_COLORS])
{
	float auxLogProbObs;
	unsigned long rowIdx = 0;
	unsigned int initStatesCount = initStates.size();
	unsigned long* outBestIdx = new unsigned long[nBeam];
	for (unsigned long i = 0; i < initStatesCount; i++)
	{
		auxLogProbObs = 0;
		rowIdx = i * N_COLORS;
		for (unsigned long j = 0; j < N_COLORS; j++)
			auxLogProbObs -= (float)SQUARE(((obs[j] - initMeansH[rowIdx + j]) / initStdsH[rowIdx + j]));
		auxLogProbObs /= 2;
		auxLogProbObs += initLogProbTermH[i];
		initStatesLogProbHCopy[i] += auxLogProbObs;
	}
	highestKValsInArray(initStatesLogProbHCopy.data(), initStatesLogProbHCopy.size(), nBeam, outBestIdx);
	for (unsigned int i = 0; i < nBeam; i++)
	{
		unsigned long state_idx = outBestIdx[i];
		(*outMostLikely)[i] = (initStates[state_idx]);
		if (i == 0)
			(*outProbsNorm)[i] = (initStatesLogProbHCopy[state_idx]);
		else
			(*outProbsNorm)[i] = (initStatesLogProbHCopy[state_idx] - outProbsNorm->front());
	}
	outProbsNorm->front() = 0; //Normalized log probs
	delete[] outBestIdx;
}

void InitStates::precomputeForHardcoreInitStates()
{
	unsigned int initStatesCount = initStates.size();
	initLogProbTermH.reserve(initStatesCount); initMeansH.reserve(N_COLORS * initStatesCount); initStdsH.reserve(N_COLORS * initStatesCount);
	unsigned long idxFor2dArrays;
	for (unsigned long i = 0; i < initStatesCount; i++)
	{
		idxFor2dArrays = i * N_COLORS;
		initMeansH.push_back(initStates[i].K[0] * MU); initMeansH.push_back(initStates[i].K[1] * MU); initMeansH.push_back(initStates[i].K[2] * MU);
		initStdsH.push_back(getStdState(initStates[i], 0)); initStdsH.push_back(getStdState(initStates[i], 1)); initStdsH.push_back(getStdState(initStates[i], 2));
		initLogProbTermH.push_back((float)-log(LOGTERM_CONST * initStdsH[idxFor2dArrays] * initStdsH[idxFor2dArrays + 1] * initStdsH[idxFor2dArrays + 2]));
	}
}

float InitStates::calcDyeLossProbLog(unsigned int Kinit[N_COLORS], unsigned int Kend[N_COLORS])
{
	float retVal = 0;
	float totalInitDyes = Kinit[0] + Kinit[1] + Kinit[2];
	float comb_factor = (float)(nChoosek(Kinit[0], Kend[0]) * nChoosek(Kinit[1], Kend[1]) * nChoosek(Kinit[2], Kend[2])); //Combinational factor of picking which dyes were not attached
	float dif = totalInitDyes - Kend[0] - Kend[1] - Kend[2];
	retVal = log(comb_factor) + dif * logL +  (totalInitDyes - dif ) * logOneMinusL;
	return retVal;
}

//Useful functions

vector<unsigned int> argsortf(const vector<float>& vf) { //Argsort for a float vector. based on https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
	vector<unsigned int> idx(vf.size());// initialize original index locations
	iota(idx.begin(), idx.end(), 0);
	stable_sort(idx.begin(), idx.end(),
		[&vf](unsigned int i1, unsigned int i2) {return vf[i1] > vf[i2]; });
	return idx;
}



unsigned int nChoosek(unsigned int n, unsigned int k) //https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
{
	if (k > n) return 0;
	if (k * 2 > n) k = n - k;
	if (k == 0) return 1;

	unsigned int result = n;
	for (unsigned int i = 2; i <= k; ++i) {
		result *= (n - i + 1);
		result /= i;
	}
	return result;
}

void highestKValsInArray(float* array, unsigned long n, unsigned int k, unsigned long* outBestIdx) //gets higher k values in a long array, where k<<n
{
	float* bestProbs = new float[k];
	int auxK = 0;
	for (unsigned long i = 0; i < k; i++)
		bestProbs[i] = -BIG_FLOAT; //Set a very low value
	for (unsigned long i = 0; i < n; i++)
	{
		if (array[i] > bestProbs[k - 1]) //If it is best than the worst best
		{
			auxK = k - 1;
			while (array[i] > bestProbs[auxK]) //Continues comparing until finding a higher prob or getting to the first ones
				if (auxK-- == 0)
					break;
			auxK++;
			for (int j = k - 2; j >= auxK; j--) { //All probs and indexes are pushed from the known ones
				bestProbs[j + 1] = bestProbs[j];
				outBestIdx[j + 1] = outBestIdx[j];
			}
			bestProbs[auxK] = array[i]; //New data is inserted in the list
			outBestIdx[auxK] = i;
		}
	}
	delete[] bestProbs;
}

