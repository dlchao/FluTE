/* class EpiModelParameters
 *
 * Dennis Chao
 * 12/2009
 */

#include <climits>
#include <cstring>
#include <string>
using namespace std;
#include "params.h"

class EpiModelParameters {
 public:
#ifdef PARALLEL
  EpiModelParameters(int, int, const char *configname);
#else
  EpiModelParameters(const char *configname);
#endif
  ~EpiModelParameters() { }

  // run parameters
  int getRunLength() { return nRunLength; }
  string getLabel() { return szLabel; }
  string getBaseName() { return szBaseName; }
  unsigned int getRandomSeed() { return seeddisp; }

  // epidemic parameters
  double getBeta() { return beta; }
  double getR0() { return R0; }
  double *getSeasonality() { return seasonality; }
  int getTriggerDelay() { return nTriggerDelay; }
  int getTriggerDay() { return nTriggerDay; }
  int getAscertainmentDelay() { return nAscertainmentDelay; }
  double getAscertainmentFraction() { return fSymptomaticAscertainment; }
  double getContactAscertainmentFraction() { return fContactAscertainment; }
  bool getSeedDaily() { return bSeedDaily; }
  int getSeedAirports() { return nSeedAirports; }
  bool getTravel() { return bTravel; }
  unsigned int getSeedInfectedTractFIPS() { return nSeedInfectedTractFIPS; }
  unsigned int getSeedInfectedCountyFIPS() { return nSeedInfectedCountyFIPS; }
  unsigned int getSeedInfectedStateFIPS() { return nSeedInfectedStateFIPS; }
  int getSeedInfectedNumber() { return nSeedInfectedNumber; }
  double *getPreexistingImmunityByAge() { return fPreexistingImmunityFraction; }
  double getPreexistingImmunityLevel() { return fPreexistingImmunityLevel; }
  double *getBaselineVESByAge() { return fBaselineVESByAge; }
  int *getSchoolOpeningDays() { return nSchoolOpeningDays; }

  // vaccination/AV parameters
  enum prevaccType getPrevaccinationStrategy() { return ePrevaccinationStrategy; }
  enum vaccType getVaccinationStrategy() { return eVaccinationStrategy; }
  unsigned char *getVaccinePriorities() { return nVaccinePriorities;}
  unsigned char *getVaccinePriorities2() { return nVaccinePriorities2;}
  int getPriorityChangeTime() { return nPriorityChangeTime;}
  double getVaccinationFraction() { return fVaccinationFraction; }
  double getResponseThreshhold() { return fResponseThreshhold; }
  unsigned int *getVaccineInitialSupply() { return nVaccineInitialSupply; }
  int getVaccineDailyLimit() { return nVaccineDailyLimit; }
  int getAVTotalLimit() { return nAVTotalLimit; }
  int getAVDailyLimit() { return nAVDailyLimit; }
  enum antiviralPolicy getAVPolicy() { return eAntiviralPolicy; }
  vaccinedatastruct *getVaccineData() { return vaccinedata; }
  double *getVaccineEfficacyByAge() { return fVaccineEfficacyByAge; }
  bool *getVaccineBoostByAge() { return bVaccineBoostByAge; }
  unsigned int *getVaccineProductionSchedule(int vnum) { return vaccineproductionschedule[vnum]; }
  double getAVEs() { return AVEs; }
  double getAVEp() { return AVEp; }
  double getAVEi() { return AVEi; }
  double getAdultEssentialFraction() {return fAdultEssentialFraction;}
  double *getPregnantFraction() {return fPregnantFraction;}
  double *getHighRiskFraction() {return fHighRiskFraction;}
 
  // NPI parameters
  int getSchoolClosurePolicy() { return schoolClosurePolicy; }
  int getSchoolClosureDays() { return nSchoolClosureDays; }
  double getIsolationCompliance() { return fIsolationCompliance; }
  double getQuarantineCompliance() { return fQuarantineCompliance; }
  double getLiberalLeaveCompliance() { return fLiberalLeaveCompliance; }

  // log files
  bool getIsLogFile() { return (nLogFileInterval>0); }
  int getLogFileInterval() { return nLogFileInterval; }
  ofstream *getTractFile() { return &tractFile; }
  ofstream *getLogFile() { return &logFile; }
  ofstream *getSummaryFile() { return &summaryFile; }
  ofstream *getIndividualsFile() { return &individualsFile; }
  bool getHasIndividualsFile() { return bIndividualsFile; }

 protected:
  bool read_config_bool(bool &result, istream &iss, const char *s);
  bool read_config_int(int &result, istream &iss, const char *s);
  bool read_config_unsigned(unsigned int &result, istream &iss, const char *s);
  bool read_config_double(double &result, istream &iss, const char *s, const double min, const double max);
  bool readConfigFile(const char *configname);
  bool createOutputFiles(void);

  int nRunLength; // total time of simulation in days
  double beta;    // transmisison parameter
  double R0;      // R_0 (used to derive beta)
  double seasonality[MAXRUNLENGTH];      // multiplier of beta
  double fPreexistingImmunityFraction[TAG];  // fraction with pre-existing (sterilizing) immunity by age group
  double fPreexistingImmunityLevel; // protection against infection for those with pre-existing immunity (like VES)
  double fBaselineVESByAge[TAG];    // default VES by age group

  // logging parameters
  int nLogFileInterval;// output to log file every n days (set to 0 for no log)
  ofstream tractFile;
  ofstream logFile;
  ofstream summaryFile;
  ofstream individualsFile;

  // run parameters
  string szLabel;       // name of this run
  int nr;	        // the run number
  int seeddisp;	        // random number generator seed displacement
  bool bIndividualsFile; // whether to output an individuals file
  string szBaseName;	// basename of data files
  string szLogFileName;      // name of log file
  string szSummaryFileName;  // name of summary file 
  string szTractFileName;    // name of tract file 
  string szIndividualFileName;    // name of individual file 
  unsigned int nSeedInfectedTractFIPS;  // FIPS of the census tract to seed infected people
  unsigned int nSeedInfectedCountyFIPS; // FIPS of the county to seed infected people
  unsigned int nSeedInfectedStateFIPS;  // FIPS of the state to seed infected people
  int nSeedInfectedNumber;     // number of infected people to seed
  bool bSeedDaily;             // seed infected people every day or just once?
  int  nSeedAirports;          // seed n/10000 infected people at major airports
  bool bTravel;                // is short-term travel allowed?
  int nSchoolOpeningDays[56];

  // intervention parameters
  int bTrigger;                     // has the trigger for response been reached
  int nTriggerDay;                  // day when reactive strategies will be deployed everywhere
  int nTriggerDelay;                // days between trigger and response
  double fResponseThreshhold;       // fraction of ever infecteds before reaction
  int nAscertainmentDelay;          // days between symptomatic and ascertainment
  double fAdultEssentialFraction;   // fraction of working-age adults who are "essential workers"
  double fPregnantFraction[TAG];    // fraction of people who are pregnant, by age
  double fHighRiskFraction[TAG];    // fraction of individuals who are at high risk for mortality/morbidity

  // pharmaceutical intervention parameters
  enum prevaccType ePrevaccinationStrategy; // actions to take before simulation starts
  enum vaccType eVaccinationStrategy;       // actions to take in response to epidemic
  unsigned char nVaccinePriorities[PRIORITY_LAST]; // vaccination priorities for all categories
  unsigned char nVaccinePriorities2[PRIORITY_LAST]; // vaccination priorities for all categories
  int nPriorityChangeTime;           // when to switch to "second" policy
  double fVaccinationFraction;      // vaccinate this proportion of the population
  unsigned int uVaccinationFraction;// unsigned int version of fVaccinationFraction
  double fContactAscertainment;         // fraction of ascertainment of non-close mixing groups
  double fSymptomaticAscertainment;    // fraction of ascertainment of symptomatic people
  enum antiviralPolicy eAntiviralPolicy; // how to distribute antivirals
  unsigned int nVaccineInitialSupply[NUMVACCINES];  // total doses of vaccines initially available
  unsigned int nAVTotalLimit;       // total courses of AV available
  unsigned int nVaccineDailyLimit;  // doses of vaccines that can be used per day (distribution capacity)
  unsigned int nAVDailyLimit;       // courses of AV that can be used per day

  vaccinedatastruct vaccinedata[NUMVACCINES];
  double fVaccineEfficacyByAge[TAG];           // relative efficacy of vaccine by age
  bool bVaccineBoostByAge[TAG+3]; // requires boost?
  unsigned int vaccineproductionschedule[NUMVACCINES][MAXRUNLENGTH];
  double AVEs;  // less likely to get infected
  double AVEp;  // less likely to become sick given infection
  double AVEi; // less likely to infect others

  // non-pharmaceutical intervention parameters
  int schoolClosurePolicy;          // 
  int nSchoolClosureDays;           // number of days to close schools (0 for no school closure)
  double fIsolationCompliance;      // probability of voluntary home isolation compliance (set to 0 for no isolation)?
  double fQuarantineCompliance;     // probability of individual compliance (set to 0 for no quarantine)
  double fLiberalLeaveCompliance;   // probability of individual compliance (set to 0 for no liberal leave)

#ifdef PARALLEL
  int size;
  int rank;
#endif
};
