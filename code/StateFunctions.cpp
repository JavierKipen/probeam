#include "StateFunctions.h"
#include <cmath>



#define ABS(x)  ( ( (x) < 0) ? -(x) : (x) )



unsigned int getCountN(State& s)
{
	unsigned int res=0;
	for (unsigned int i = 0; i < N_COLORS; i++)
		res += s.N[i];
	return res;
}

unsigned int getCountK(State& s)
{
	unsigned int res = 0;
	for (unsigned int i = 0; i < N_COLORS; i++)
		res += s.K[i];
	return res;
}

bool finishedSequencing(State& s)
{
	bool res = true;
	for (unsigned int i = 0; i < N_COLORS; i++)
		res &= (s.K[i]==0);
	return res;
}

bool isEqual(State& s1, State& s2)
{
	bool retVal = false;
	bool resN, resK, resR;
	resK = true;
	for (unsigned int i = 0; i < N_COLORS; i++)
		resK &= (s1.K[i] == s2.K[i]);
	if (resK)
	{
		resN = true;
		for (unsigned int i = 0; i < N_COLORS; i++)
			resN &= (s1.N[i] == s2.N[i]);
		if (resN && (s1.RCharCount == s2.RCharCount))
		{
			resR = true;
			for (unsigned int i = 0; i < s1.RCharCount; i++)
				resR &= (s1.R[i] == s2.R[i]);
			if (resR)
				retVal = true;
		}
	}
	return retVal;
}
bool isNEqual(State& s1, unsigned int Np[N_COLORS])
{
	bool retVal = true;
	for (unsigned int i = 0; i < N_COLORS; i++)
		retVal &= (s1.N[i] == Np[i]);
	return retVal;
}

float pDyeAttached(State& s, unsigned int i)
{
	return (s.N[i]==0) ? 0: (float)s.K[i] / (float)s.N[i];
}
float pDyeAttached(unsigned int N[N_COLORS], unsigned int K[N_COLORS], unsigned int i)
{
	return (N[i] == 0) ? 0 : (float)K[i] / (float)N[i];
}

float logpObs(State& s, float obs[N_COLORS])
{
	float means_loc[N_COLORS];
	float stds_loc[N_COLORS];
	float retVal = 0; //Logprob of 1;

	for (unsigned long j = 0; j < N_COLORS; j++) //Obtains values of means and stds with the given K
	{
		means_loc[j] = getMeanState(s, j); 
		stds_loc[j] = getStdState(s, j);
		retVal -= (float)pow((obs[j] - means_loc[j]) / stds_loc[j], 2);
	}
	retVal /= 2; //Division of 2 because of the exponent
	retVal -= (float)log(LOGTERM_CONST * stds_loc[0] * stds_loc[1] * stds_loc[2]); //constant term
	return retVal;
}

float estExpTerm(unsigned int K[N_COLORS], float obs[N_COLORS])
{
	float retVal = 0;
	for (unsigned long j = 0; j < N_COLORS; j++) //Obtains values of means and stds with the given K
	{
		if (K[j] == 0)
			retVal += SQUARE(obs[j] - MU * K[j]) / STD_B2;
		else
			retVal += SQUARE(obs[j] - MU * K[j]) / (K[j]* STD2);
	}
	return retVal;
}
State* decK(State* s, unsigned int i)
{
	s->K[i] -= 1;
	return s;
}

State* decN(State* s, unsigned int i)
{
	s->N[i] -= 1;
	return s;
}

State* detach(State* s)
{
	unsigned int zeros[N_COLORS] = { 0,0,0 }; //This could be written in a non hardcoded way
	setK(s, zeros);
	return s;
}

State* setK(State* s, unsigned int Kp[N_COLORS])
{
	for (unsigned int i = 0; i < N_COLORS; i++)
		s->K[i] = Kp[i];
	return s;
}

State* setN(State* s, unsigned int Np[N_COLORS])
{
	for (unsigned int i = 0; i < N_COLORS; i++)
		s->N[i] = Np[i];
	return s;
}

State* copyState(State* dest, State& orig)
{
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		dest->N[i] = orig.N[i];
		dest->K[i] = orig.K[i];
	}
	for (unsigned int i = 0; i < orig.RCharCount; i++)
	{
		dest->R[i] = orig.R[i];
	}
	for (unsigned int i = 0; i < orig.dyeSeqsIdxsCount; i++)
	{
		dest->dyeSeqsIdxs[i] = orig.dyeSeqsIdxs[i];
	}
	dest->dyeSeqsIdxsCount = orig.dyeSeqsIdxsCount;
	dest->RCharCount = orig.RCharCount;
	return dest;
}

State* remChar(State* s, char ch)
{
	s->R[s->RCharCount++] = ch;
	return s;
}
State* setNewDyeSeqsIdxs(State* s, unsigned int* dyeSeqsIdxs, unsigned int newDyeSeqsIdxsCount)	//Copies new dye seqs idxs
{
	for (unsigned int i = 0; i < newDyeSeqsIdxsCount; i++)
		s->dyeSeqsIdxs[i] = dyeSeqsIdxs[i];
	s->dyeSeqsIdxsCount = newDyeSeqsIdxsCount;
	return s;
}


float getStdState(State& S, unsigned int stdIdx)
{
	float std = STD_B;
	std = sqrt((S.K[stdIdx] * STD2) + STD_B2);
	return std;
}
float getMeanState(State& S, unsigned int meanIdx)
{
	float mean = 0;
	mean = S.K[meanIdx] * MU;
	return mean;
}


float probObsSingleColor(float x, unsigned int k) //https://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
{
	float std = sqrt((k * STD2) + STD_B2);
	float mu = k * MU;
	static const float inv_sqrt_2pi = 0.3989422804014327;
	float aux = (x - mu) / std;

	return inv_sqrt_2pi / std * exp(-0.5f * aux * aux);
}

float logpObs(unsigned int K[N_COLORS], float obs[N_COLORS])
{
	float means_loc[N_COLORS];
	float stds_loc[N_COLORS];
	float retVal = 0; //Logprob of 1;
	float auxvar;

	for (unsigned long j = 0; j < N_COLORS; j++) //Obtains values of means and stds with the given K
	{
		means_loc[j] = K[j] * MU;
		stds_loc[j] = sqrt((K[j] * STD2) + STD_B2);
		auxvar = (obs[j] - means_loc[j]) / stds_loc[j];
		retVal -= (float) auxvar * auxvar;
	}
	retVal /= 2; //Division of 2 because of the exponent
	retVal -= (float)log(LOGTERM_CONST * stds_loc[0] * stds_loc[1] * stds_loc[2]); //constant term
	return retVal;
}

bool isNEqual(unsigned int N1[N_COLORS], unsigned int N2[N_COLORS])
{
	bool retVal = true;
	for (unsigned int i = 0; i < N_COLORS; i++)
		retVal &= (N1[i] == N2[i]);
	return retVal;
}

bool isNZero(unsigned int N[N_COLORS])
{
	unsigned int Z[N_COLORS] = { 0 };
	return isNEqual(N, Z);
}

unsigned int KDist(unsigned int K1[N_COLORS], unsigned int K2[N_COLORS])
{
	unsigned int retVal = 0;
	for (unsigned int i = 0; i < N_COLORS; i++)
		retVal += ABS( (int)K1[i] - (int)K2[i]);
	return retVal;
}



unsigned int getBestKi(float obs)
{
	unsigned int auxKi = (int(obs / MU) < 0) ? 0 : int(obs / MU); //Lowest possible value is floor or 0.
	auxKi = (probObsSingleColor(obs, auxKi) > probObsSingleColor(obs, auxKi + 1)) ? auxKi : auxKi + 1; //Keeps the most likely Ki
	return auxKi;
}

bool isK1leqK2elem(unsigned int K1[N_COLORS], unsigned int K2[N_COLORS])
{
	bool retVal=true;
	for (unsigned int i = 0; i < N_COLORS; i++)
		retVal &= (K1[i] <= K2[i]);
	return retVal;
}

StateRed* copyState(StateRed* dest, State& orig)
{
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		dest->N[i] = orig.N[i];
		dest->K[i] = orig.K[i];
	}
	for (unsigned int i = 0; i < orig.RCharCount; i++)
	{
		dest->R[i] = orig.R[i];
	}
	dest->RCharCount = orig.RCharCount;
	return dest;
}
StateRed* setK(StateRed* s, unsigned int Kp[N_COLORS])
{
	for (unsigned int i = 0; i < N_COLORS; i++)
		s->K[i] = Kp[i];
	return s;
}

StateRed* remChar(StateRed* s, char ch)
{
	s->R[s->RCharCount++] = ch;
	return s;
}
StateRed* decN(StateRed* s, unsigned int i)	//Decreases K at index i(MODIFIES STATE)
{
	s->N[i] -= 1;
	return s;
}

StateRed* omitRemoval(StateRed* s)
{
	s->RCharCount--;
	return s;
}
StateRed* incN(StateRed* s, unsigned int i)
{
	s->N[i] += 1;
	return s;
}
StateRed* decK(StateRed* s, unsigned int i)
{
	s->K[i] -= 1;
	return s;
}
StateRed* incK(StateRed* s, unsigned int i)
{
	s->K[i] += 1;
	return s;
}
unsigned int* incVect(unsigned int* V, unsigned int i)
{
	V[i] += 1;
	return V;
}
unsigned int* decVect(unsigned int* V, unsigned int i)
{
	V[i] -= 1;
	return V;
}
bool isEqual(StateRed& s1, StateRed& s2)
{
	bool retVal = false;
	bool resN, resK, resR;
	resK = true;
	for (unsigned int i = 0; i < N_COLORS; i++)
		resK &= (s1.K[i] == s2.K[i]);
	if (resK)
	{
		resN = true;
		for (unsigned int i = 0; i < N_COLORS; i++)
			resN &= (s1.N[i] == s2.N[i]);
		if (resN && (s1.RCharCount == s2.RCharCount))
		{
			resR = true;
			for (unsigned int i = 0; i < s1.RCharCount; i++)
				resR &= (s1.R[i] == s2.R[i]);
			if (resR)
				retVal = true;
		}
	}
	return retVal;
}

void copyState(State* dest, StateRed& orig) // Copies without dye sequences
{
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		dest->N[i] = orig.N[i];
		dest->K[i] = orig.K[i];
	}
	for (unsigned int i = 0; i < orig.RCharCount; i++)
	{
		dest->R[i] = orig.R[i];
	}
	dest->RCharCount = orig.RCharCount;
}

float approxDistZScore(unsigned int Kp[N_COLORS], float obs[N_COLORS])
{
	float retVal = 0;
	float aux;
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		if (Kp[i] == 0)
			aux = ABS(obs[i] / STD_B);
		else
			aux = ABS((obs[i] - MU * Kp[i]) / (sqrt(Kp[i])*STD));
#ifdef PRUN_WORST_Z_SCORE
		if (aux > retVal)
			retVal = aux;
#else
		retVal += aux;
#endif // 

		
	}
	return retVal;
}