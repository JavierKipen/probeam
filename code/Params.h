#pragma once

#define N_COLORS 3
#define N_ED_CYC 9
#define N_FEATURES_PER_COL (N_ED_CYC+1)
#define N_FEATURES_TOTAL (N_FEATURES_PER_COL*3)

#define E 0.06f
#define M 0.07f
#define L 0.05f
#define D 0.05f
#define MU 10000.0f
#define STD 1600.0f
#define STD2 (STD*STD)
#define STD_B 66.0f
#define STD_B2 (STD_B*STD_B)

#define SQUARE(i) ( (i) * (i) )

#define ENUMBER (2.71828)
#define LOGTERM_CONST 15.7496099457 //(2*PI)^(3/2) precomputed to make it faster

#define DEF_NBEAM 30
#define ZSCOREDEFAULT 5
#define PRUN_WORST_Z_SCORE //(FOR PRUNING)If not commented looks at the worst z score, if commented looks at the sum of z scores 

//#define ESTIMATE_IFED_TIMES //When defined, keeps track of the timings for the calculations of IFEDs.