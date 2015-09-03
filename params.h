/*
 * params.h
 * some constants and data structures for FluTE
 */

#ifndef __PARAMS_H
#define __PARAMS_H

enum {
  AG0,  // 0-4 years old
  AG1,  // 5-18 years old
  AG2,  // 19-29 years old
  AG3,  // 30-64 years old
  AG4,	// 65 and older
  TAG	// total age groups
};	// age groups

extern const double fAG0InfantFraction; // fraction of those in age group 0 who are <6 months

// contact probabilities
extern const double cpw;       // workplace 
extern const double cps[10];   // schools/daycare
extern const double cpcm[TAG]; // community
extern const double cpnh[TAG]; // neighborhood
extern const double cpfc[TAG]; // family from child
extern const double cpfa[TAG]; // family from adults
extern const double cphcc[TAG];// household cluster from child
extern const double cphca[TAG];// household cluster from adult

// withdraw probabilities
#define WITHDRAWDAYS 3
extern const double withdrawprob[3][WITHDRAWDAYS];

// self isolation probabilities by age
extern const double isolationprob[TAG];

extern const int nQuarantineLength; // length of quarantine in days
extern const int nAntiviralCourseSize; // number of pills in one antiviral course (1 tablet/day for prophylaxis and 2 tablets/day for treatment)
extern const double fStopAntiviralTwoPills; // probability that individuals taking antivirals stop after exactly two pills

// cdf of duration for incubation period in days
extern const double incubationcdf[3];
//extern const unsigned int incubationcdf32[3];

// viral load trajectories
// number of subjects
#define VLOADNSUB	6
// number of days with detectable viral load
#define VLOADNDAY	6
extern const double basevload[VLOADNSUB][VLOADNDAY];

// vaccine efficacy over time
// 28 days for VE buildup
#define VACCEFFLENGTH 28

// maximum length of a simulation in days (change this if you want longer runs)
#define MAXRUNLENGTH 365

// travel data
extern const double travel_pr[TAG]; // age-specific travel probability per day
extern const double travel_length_cdf[12]; // cdf of length of trip in days (nights away from home)

#define NUMBER_HUBS 15
extern const unsigned int FIPS_hubs[];
extern const unsigned int Size_hubs[];
// 15 busiest airports from US International Air Passenger and Freight Statistics, June 2008
#endif
