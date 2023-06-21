#pragma once
#include "Params.h"
#include "State.h"

using namespace std;

class InfoForEdmanDegradation {
public:
	InfoForEdmanDegradation();
	void clear();
	bool dyeCanBeRemoved[N_COLORS];
	bool nonLumCanBeRemoved;
	float pRemDye[N_COLORS];
	float pRemNonLum;
	unsigned int dyeSeqsIdxs[N_COLORS][N_MAX_DYESEQS_IN_STATE];	//Next ids for all posibilities of removal of color
	unsigned int dyeSeqsIdxsCount[N_COLORS];					//Count of ids for every possible removal of color
	unsigned int dyeSeqsIdxsDot[N_MAX_DYESEQS_IN_STATE]; // Next ids for possibilities of removal of .
	unsigned int dyeSeqsIdxsCountDot; // Next ids for possibilities of removal of .
};

typedef InfoForEdmanDegradation IFED;

