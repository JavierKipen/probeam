#include "InfoForEdmanDegradation.h"


InfoForEdmanDegradation::InfoForEdmanDegradation()
{
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		dyeCanBeRemoved[i] = false;
		pRemDye[i] = 0;
		dyeSeqsIdxsCount[i] = 0;
	}
	dyeSeqsIdxsCountDot = 0;
	nonLumCanBeRemoved = false;
	pRemNonLum = 0;
}

void InfoForEdmanDegradation::clear()
{
	for (unsigned int i = 0; i < N_COLORS; i++)
	{
		dyeCanBeRemoved[i] = false;
		pRemDye[i] = 0;
		dyeSeqsIdxsCount[i] = 0;
	}
	dyeSeqsIdxsCountDot = 0;
	nonLumCanBeRemoved = false;
	pRemNonLum = 0;
}