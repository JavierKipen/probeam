#pragma once
#include "Params.h"

#define N_MAX_DYESEQS_IN_STATE 7500 // Should be obtained in the worst case (whole proteome)

typedef struct {
	unsigned int N[N_COLORS];	//Number of ideal fluorophores remaining for each color
	unsigned int K[N_COLORS];	//Fluorophores remaining
	char R[N_ED_CYC];			//Removed sequence chars
	unsigned int RCharCount;	//Count of the chars
	unsigned int dyeSeqsIdxs[N_MAX_DYESEQS_IN_STATE];	//Idxs of dye sequences that belong to this state
	unsigned int dyeSeqsIdxsCount; //Count of the dye sequences
} State;
