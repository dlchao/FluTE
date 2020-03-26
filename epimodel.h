/* class EpiModel
 *
 * This class models influenza transmission among people in a
 * synthetic population as described in the paper "FluTE, a publicly 
 * available stochastic influenza epidemic simulation model," 
 * available at http://www.ploscompbiol.org/doi/pcbi.1000656
 *
 * Dennis Chao and Shufu Xu
 * 12/2009
 */

#include <iterator>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <assert.h>
#include "params.h"
extern "C" {
  #include "dSFMT.h"   // for random number generation
}

using namespace std;

class EpiModelParameters;

enum {
	HELLO,
	DATA
};	// MPI tags

// pre-vaccination and priming options
enum prevaccType {
  NOPREVACC=0,
  PREVACCINATE,  // pre-vaccinate the population
  PRIMEBOOSTRANDOM,// prime the population with one dose, boost random people later
  PRIMEBOOSTSAME   // prime the population with one dose, boost same people later
};

// dynamic vaccination options
enum vaccType {
  NOVACC=0,
  RINGVACTRACT,
  RINGVACCOUNTY,
  MASSVAC
};

enum vaccinePriorityPolicy {
  PRIORITY_ESSENTIAL=0,    // essential workforce
  PRIORITY_PREGNANT=1,     // pregnant women
  PRIORITY_INFANTFAMILY=2, // family members of infants
  PRIORITY_HR0=3,          // high risk preschoolers
  PRIORITY_HR1=4,          // high risk school-age children
  PRIORITY_HR2=5,          // high risk young adults
  PRIORITY_HR3=6,          // high risk older adults
  PRIORITY_HR4=7,          // high risk elderly
  PRIORITY_0=8,            // low risk preschoolers
  PRIORITY_1=9,            // low risk school-age children
  PRIORITY_2=10,           // low risk young adults
  PRIORITY_3=11,           // low risk older adults
  PRIORITY_4=12,           // low risk elderly
  PRIORITY_LAST=13
};

// antiviral distribution policy options
enum antiviralPolicy {
  NOAV=0,
  TREATMENTONLY=1,
  HHTAP=2,
  HHTAP100=3,
  FULLTAP=4
};

enum {
  NoError,
  MemError,
  FileError,
  RunfileError,
  LogfileError,
  TractfileError,
  SummaryError,
  SyncError,
  MiscError
};	// errors

enum {
  SUSCEPTIBLE = 0x1u,
  INFECTED    = 0x2u,
  SYMPTOMATIC = 0x4u,
  VACCINATED  = 0x8u,  // at least one dose of vaccine
  BOOSTED     = 0x10u, // received both doses
  ANTIVIRAL   = 0x20u, // at least one dose of antiviral
  NEEDSBOOST  = 0x40u  // needs a 2nd dose of a 1-dose vaccine
};      // status bits for persons (status)

enum {
  WILLBESYMPTOMATIC = 0x40u, // will become symptomatic if infected
  WILLBEASCERTAINED = 0x80u // will become ascertained if infected
};      // status bits for persons (ibits). Note: the first few bits of ibits are used to keep track of incubation and withdrawal timers

enum {
  //  VACCINE2       = 0x1u, // first 4 bits for vaccine ID
  WITHDRAWN      = 0x20u,
  QUARANTINED    = 0x40u,
  AVPROPHYLAXIS  = 0x80u  // given antivirals for prophylaxis (not symptomatic)
  //  ONLIBERALLEAVE = 0x8u (not needed)
  //  ISOLATED       = 0x4u, (not needed)
};      // more status bits for persons (vbits)

enum vrestrictionbits {
  VESSENTIAL = 0x01u,   // essential worker
  VPREGNANT  = 0x02u,   // pregnant 
  VHIGHRISK  = 0x04u,   // high risk
  VINFANT    = 0x08u    // <6 months old
}; // vaccination restriction/priority status

enum tractstatusbits {
  TRACTAV1 = 0x1u,         // AV policy used in tract (low order bit)
  TRACTAV2 = 0x2u,         // AV policy used in tract (low order bit)
  LIBERALLEAVE = 0x4u,     // liberal leave policy for workers
  TRACTVACCINATED = 0x8u, // tract-wide vaccination
  TRACTQUARANTINE = 0x10u  // household quarantine active in this tract
  //  TRAVELRESTRICTED = 0x20u, // reduced long-distance travel
};      // status bits for census tracts

enum schoolclosurebits {
  DAYCARECLOSED = 0x1u,
  ELEMENTARYCLOSED = 0x2u,
  MIDDLECLOSED = 0x4u,
  HIGHCLOSED = 0x8u
};

enum {
  FROMFAMILY       = 1,
  FROMCOMMUNITY    = 2,
  FROMNEIGHBORHOOD = 3,
  FROMSCHOOL       = 4,
  FROMWORK         = 5,
  FROMFAMILYNIGHT       = 10,
  FROMCOMMUNITYNIGHT    = 11,
  FROMNEIGHBORHOODNIGHT = 12,
  FROMHHCLUSTERNIGHT    = 13,
  FROMSEED         = 20,
  FROMAIRPORT      = 21
}; // source of infection

// number of possible vaccines
#define NUMVACCINES 16

// information required for each vaccine
struct vaccinedatastruct {
  double VEs,VEi,VEp;
  int nNumDoses;     // how many doses required? (1 or 2)
  int nBoostDay;     // how many days between prime and boost? (minimum days)
  double fInfants;   // bad for this fraction of infants
  double fAge[TAG];  // bad for this fraction of age group
  bool bPregnant;    // bad for pregnant women?
  double vacceff[VACCEFFLENGTH+1]; // vaccine efficacy over time after the first shot. If a boost is required, then the vaccine efficacy will stop at nBoostDay until the boost is received.
};

struct Person {
  unsigned int id;	// agent id
  unsigned char age;	// age group
  unsigned char status;	// infection status
  unsigned char nWhichVload;    // which viral load trajectory (vload) to use?
  char iday;            // infected days
  unsigned char ibits;  // state info for infection
  unsigned char vbits;  // state info for vaccine
  unsigned char vday;   // vaccination timer (number of days elapsed)
  char nAVTimer;        // antiviral timer (number of tablets left)
  unsigned char nQuarantineTimer; // quarantine timer (number of days left)
  char nTravelTimer;	// travel time left
  char nDayNeighborhood; // ID of work neighborhood
#ifdef PARALLEL
  int nHomeRank;        // node of home community
  int nWorkRank;        // node of work community
#endif
  unsigned int nDayTract;        // ID of daytime (work, school) tract
  unsigned int nDayComm;         // ID of daytime (work, school) community
  int nWorkplace;	// work or school group
  double fBaselineVES;    // baseline VES before vaccination
  unsigned int nHomeComm;        // ID of home community

  char sourcetype;      // in which setting infection took place
  unsigned int sourceid;// keeping track of the source of infection
#ifdef PARALLEL
  int sourcerank;       // keeping track of the source of infection
#endif

  int family;	        // family ID
  unsigned char nFamilySize; // family size
  char nHomeNeighborhood;    // neighborhood (0-3)
  unsigned char householdcluster;// household cluster
  unsigned char nVaccinePriority; // vaccine priority group (0 for no vaccine, 1 for highest priority, 2 for next-highest priority, etc). Also set to 0 if the person does not want vaccine.
  int nInfectedTime;	// keeping track of the time of infection
  unsigned int nVaccineRestrictionBits; // for keeping track of vaccine restriction categories
  bool bVaccineEligible[NUMVACCINES];

  // temporary data
  double pri; // infectiousness multiplier
  double prs; // susceptibility multiplier
  bool bWantVac; // Should be vaccinated now (if available)
  bool bWantAV;  // Should get antiviral now (if available)
  friend ostream& operator<<(ostream& os, Person& p) {
    os << p.id << " " << (int)p.age << " " << p.nHomeComm << " " << p.nDayComm <<  " " << p.family << " " << (int)p.nWorkplace << " " << (int)p.nHomeNeighborhood << " " << (int)p.nFamilySize << " " << p.sourceid << " " << p.nInfectedTime << " " << (int)p.sourcetype << " ";
    if (isSusceptible(p)) {
      os << "s";
    }
    if (isInfected(p)) {
      os << "i";
    }
    if (isSymptomatic(p)) {
      os << "S";
    }
    return os;
  }

  // person age accessor functions
  friend inline bool isChild(const Person &p) { return (p.age<=1); } // is a child?
  friend inline bool isWorkingAge(const Person &p) { return (p.age==2 || p.age==3); } // old enough to work?
  
  // person.status accessor functions (vaccination and illness state)
  friend inline bool isSusceptible(const Person &p) { return p.status&SUSCEPTIBLE; }
  friend inline void setSusceptible(Person &p) {p.status|=SUSCEPTIBLE;} 
  friend inline void clearSusceptible(Person &p) {p.status&=~SUSCEPTIBLE;}
  friend inline bool isInfected(const Person &p) { return p.status&INFECTED; }
  friend inline bool isInfectious(const Person &p) { return (p.status&INFECTED && p.iday>=0); }
  friend inline bool isVaccinated(const Person &p) { return p.status&VACCINATED; }
  friend inline bool isBoosted(const Person &p) { return p.status&BOOSTED; }
  friend inline bool isAntiviral(const Person &p) { return p.status&ANTIVIRAL; }
  friend inline bool isSymptomatic(const Person &p) { return p.status&SYMPTOMATIC; }
  friend inline void setSymptomatic(Person &p) {p.status|=SYMPTOMATIC;}
  friend inline void setInfected(Person &p) {p.status|=INFECTED;}
  friend inline void setVaccinated(Person &p) {p.status|=VACCINATED;}
  friend inline void setBoosted(Person &p) {p.status|=BOOSTED;}
  friend inline void setAntiviral(Person &p) {p.status|=ANTIVIRAL;}
  friend inline void clearAntiviral(Person &p) {p.status&=(~ANTIVIRAL);}
  friend inline bool needsBoost(const Person &p) { return p.status&NEEDSBOOST; }
  friend inline void setNeedsBoost(Person &p) {p.status|=NEEDSBOOST;}
  friend inline void clearNeedsBoost(Person &p) {p.status|=NEEDSBOOST;}

  // person.ibits accessor functions (response to infection)
  friend inline bool getWillBeSymptomatic(const Person &p) { return p.ibits&WILLBESYMPTOMATIC; }
  friend inline bool getWillBeAscertained(const Person &p) { return p.ibits&WILLBEASCERTAINED; }
  friend inline int getIncubationDays(const Person &p) { return p.ibits&0x7u; }
  friend inline unsigned int getWithdrawDays(const Person &p) { return (p.ibits>>3)&0x7u; } // if 0, then will not withdraw.  otherwise, withdraw this may days after infection
  friend inline void setWillBeSymptomatic(Person &p) {p.ibits|=WILLBESYMPTOMATIC; }
  friend inline void setWillBeAscertained(Person &p) {p.ibits|=WILLBEASCERTAINED; }
  friend inline void setIncubationDays(Person &p, unsigned int x) {p.ibits=(p.ibits&0xF8u)|x; }
  friend inline void setWithdrawDays(Person &p, unsigned int x) {p.ibits=(p.ibits&0xC7u)|(x<<3); }

  // person.vbits accessor functions (miscellaneous state)
  friend inline void setWhichVaccine(Person &p, unsigned char nWhich) { assert(nWhich<16); p.vbits|=nWhich; } // set which vaccine did this person got (0-15, where 0 is the default vaccine)
  friend inline unsigned char whichVaccine(const Person &p) { return p.vbits&0x0F; }	// which vaccine did this person get?
  friend inline bool isWithdrawn(const Person &p) { return p.vbits&WITHDRAWN; }
  friend inline void setWithdrawn(Person &p) {p.vbits|=WITHDRAWN;}
  friend inline void clearWithdrawn(Person &p) {p.vbits&=(~WITHDRAWN);}
  friend inline bool isQuarantined(const Person &p) { return p.vbits&QUARANTINED; }
  friend inline void setQuarantined(Person &p) {p.vbits|=QUARANTINED;}
  friend inline void clearQuarantined(Person &p) {p.vbits&=(~QUARANTINED);}
  friend inline bool isAVProphylaxis(const Person &p) { return p.vbits&AVPROPHYLAXIS; }
  friend inline void setAVProphylaxis(Person &p) {p.vbits|=AVPROPHYLAXIS;}
  friend inline void clearAVProphylaxis(Person &p) {p.vbits&=~AVPROPHYLAXIS;}

  // vaccine restriction status accessor functions
  friend inline void setEssential(Person &p) { p.nVaccineRestrictionBits|=VESSENTIAL; }
  friend inline bool isEssential(const Person &p) { return p.nVaccineRestrictionBits&VESSENTIAL; }
  friend inline void setPregnant(Person &p) { p.nVaccineRestrictionBits|=VPREGNANT; }
  friend inline bool isPregnant(const Person &p) { return p.nVaccineRestrictionBits&VPREGNANT; }
  friend inline void setHighRisk(Person &p) { p.nVaccineRestrictionBits|=VHIGHRISK; }
  friend inline bool isHighRisk(const Person &p) { return p.nVaccineRestrictionBits&VHIGHRISK; }
  friend inline void setInfant(Person &p) { p.nVaccineRestrictionBits|=VINFANT; }
  friend inline bool isInfant(const Person &p) { return p.nVaccineRestrictionBits&VINFANT; }

  friend inline unsigned char getVaccinePriority(Person &p) { return p.nVaccinePriority; }
};

#ifdef PARALLEL
typedef struct PStub { // The first few members of struct Person
  unsigned int id;	// agent id
  unsigned char age;	// age group
  unsigned char status;	// infection status
  unsigned char nWhichVload; //
  char iday;            // infected days
  unsigned char ibits;  // state info for infection
  unsigned char vbits;  // state info for vaccine
  unsigned char vday;   // vaccination timer (number of days elapsed)
  char nAVTimer;        // antiviral timer (number of tablets left)
  unsigned char nQuarantineTimer; // quarantine timer (number of days left)
  char nTravelTimer;	// travel time left
  char nDayNeighborhood; // ID of work neighborhood
  int nHomeRank;        // node of home community
  int nWorkRank;        // node of work community
  unsigned int nDayTract;        // tract of work community
  unsigned int nDayComm;         // ID of work community
  int nWorkplace;	// work group
  double fBaselineVES; // baseline VES before vaccination
} PersonStub; // for setting up people who work on different nodes

typedef struct PStatusStub {
  unsigned int id;	// agent id
  int nHomeRank;        // node of home community
  unsigned int nDayComm;// ID of work community
  unsigned char status;	// infection status
  char iday;            // infected days
  unsigned char ibits;  // state info for infection
  unsigned char vbits;  // state info for vaccine
  unsigned char vday;   // vaccination timer (number of days elapsed)
  char nTravelTimer;	// travel time left
  char sourcetype;
  unsigned int sourceid;
} PersonStatusStub; // for updating the status of workers across nodes
#endif

struct Community {
  const static int TARGETCOMMUNITYSIZE=2000;
  const static int FAMILIESPERCLUSTER = 4;
  const static int WORKGROUPSIZE = 20;

  vector< unsigned int > workers;   // vector of non-resident worker ids
  // visitors is a copy of the original individuals, so we can assign it
  // new communities, families, workgroups, etc.
  // Just keep the ID and the home rank intact.
  list <Person> visitors;           // list of short-term travelers to this community
#ifdef PARALLEL
  vector <Person> immigrantworkers; // vector of workers from other nodes
#endif
  unsigned int id;	 // unique community ID
  int nTractID;          // ID of its tract
  int nNumResidents;	 // number of residents (population size)
  unsigned int nFirstPerson, nLastPerson; // ID of first and last person(+1) in community
  int nNumWorkers;	 // number of workers (includes immigrants)
  int nNumNonWorkers;	 // number of resident non-workers
  int nNumWorkersLeaving;// number of residents who work outside the community
  int nNumWorkGroups;	 // number of work groups
  int nNumAge[TAG];	 // number of residents in different age groups
  int ninf[TAG];	 // number of residents currently infected
  int nsym[TAG];	 // number of residents currently symptomatic
  int nEverInfected[TAG];// number of residents ever infected
  int nEverSymptomatic[TAG]; // number of residents ever symptomatic
  int nEverAscertained[TAG]; // number of residents ever ascertained
  double cpcm[5];        // community-specific community contact rates
  double cpnh[5];        // community-specific neighborhood contact rates
  double daycpcm[5];     // community-specific daytime community contact rates
  double daycpnh[5];     // community-specific daytime neighborhood contact rates
  double cps[10];         // community-specific school contact rates
};

struct Tract {
  unsigned int id;          // unique tract identifier
  unsigned int fips_tract;  // census tract FIPS
  unsigned int fips_county; // county FIPS
  unsigned int fips_state;  // state FIPS
  unsigned int censuspopulation;  // population (from data file)
  unsigned int employed;    // number of employed people (from data file)
  unsigned int workforce;   // number of working age people (from data file)
  double fEmploymentProb;   // employed/workforce
  double lat;	            // latitude
  double lon;               // longitude
  unsigned int status;      // status (infected, etc)
  bool bSchoolClosed[9];    // school closure status
  unsigned int nFirstCommunity;  // ID of the first community in this tract
  unsigned int nLastCommunity;   // ID of the last community (+1) in this tract
  int nNumResidents;        // number of residents (population size)
  int nSchoolClosureTimer;
  
  friend istream& operator>>(istream& is, Tract& t) {
    char p; // because it is a comma separated file
    is >> t.fips_state >> p >> t.fips_county >> p  >> t.fips_tract >> p  >> t.censuspopulation >> p   >> t.lat >> p >> t.lon;
    return is;
  }
  friend inline enum antiviralPolicy getAVPolicy(const Tract &t) { return (enum antiviralPolicy)(t.status&0x03); }
  friend inline void setAVPolicy(Tract &t, enum antiviralPolicy p) { t.status|=p;} 
  friend inline bool isQuarantine(const Tract &t) { return t.status&TRACTQUARANTINE; }
  friend inline bool isVaccinated(const Tract &t) { return t.status&TRACTVACCINATED; } 
  friend inline bool isSchoolClosed(const Tract &t, int which) { return t.bSchoolClosed[which]; }
  friend inline void setSchoolClosed(Tract &t, int which) { t.bSchoolClosed[which]=true; }
  friend inline void setSchoolOpen(Tract &t, int which) { t.bSchoolClosed[which]=false; }
  friend inline void setQuarantine(Tract &t) {t.status|=TRACTQUARANTINE;}
  friend inline void setVaccinated(Tract &t) {t.status|=TRACTVACCINATED;}
};

class EpiModel {
 public:
#ifdef PARALLEL
  EpiModel(int, int, EpiModelParameters &params);
#else
  EpiModel(EpiModelParameters &params);
#endif
  virtual ~EpiModel() { 
    if (tractToFIPStract)
      delete(tractToFIPStract);
    if (tractToFIPScounty)
      delete(tractToFIPScounty);
    if (tractToFIPSstate)
      delete(tractToFIPSstate);
#ifdef PARALLEL
    if (tractToCPU)
      delete(tractToCPU);
    if (numTractsPerNode)
      delete(numTractsPerNode);
    if (numPeoplePerNode)
      delete(numPeoplePerNode);
#endif
  }
  virtual void prerun(void);
  virtual void run(void);

 protected:
  virtual void seedinfected(void);
  void read_tracts(void);
  bool read_workflow_file(string s, unsigned int *flow, map<int,int> *statetracts);
  void read_workflow(void);
  void create_families(Community& comm, int nTargetSize);
  void create_person(int nAgeGroup, int nFamilySize, int nFamily, int nHouseholdCluster, int nNeighborhood, int nSchoolgroup, Community& comm);
  void initPopulation(void);
  bool isAscertained(const Person &p);
  bool infect(Person& p, const Person& source, double baseprob, int sourcetype);
  virtual void infect(Person& p);
  void vaccinate(Person& p);
  void vaccinate(Tract& t);
  void TAP(Person& p);
  void dayinfectsusceptibles(const Person &infected, Community &comm);
  void nightinfectsusceptibles(const Person &infected, Community &comm);
  void day(void);
  void night(void);
  void travel_start(void);
  void travel_end(void);
  void response(void);
  virtual void log(void);
  virtual void summary(void);
  virtual void outputIndividuals(void);
  bool isEligible(const Person &p, int nVacNum);

#ifdef PARALLEL
  void sync(void);
#endif

 protected:
  int nRunLength;	        // total time of simulation in days
  int nTimer;	                // number of 12-hour steps progressed
  unsigned int nNumPerson;      // number of people on this node
  unsigned int nNumFamilies;    // number of families on this node
  unsigned int nNumCommunities; // number of communities on this node
  unsigned int nFirstTract;     // id of the first tract on this node
  unsigned int nLastTract;      // id of the last tract (+1) on this node
  unsigned int nNumTractsTotal; // number of tracts across all CPUs
  double beta; // transmisison parameter
  double R0;   // R_0 (used to derive beta)
  double seasonality[MAXRUNLENGTH];      // multiplier of beta

  vector < Person >    pvec;	   // vector of residents
  vector < Tract >     tractvec;   // the census tracts
  vector < Community > commvec;    // the communities
  unsigned int *tractToFIPStract;  // census tract FIPS
  unsigned int *tractToFIPScounty; // county FIPS
  unsigned int *tractToFIPSstate;  // state FIPS
  double withdrawcdf[3][WITHDRAWDAYS];  // withdraw to home probabilities cdf

  dsfmt_t dsfmt;                   // random number generator

  // run parameters
  string szLabel;       // name of this run
  int seeddisp;	        // random number generator seed displacement
  bool bTractFile;      // whether to output a tract file
  bool bIndividualsFile;// whether to output an individuals file
  string szBaseName;	// basename of data file

  // epidemic parameters
  unsigned int nSeedInfectedTractFIPS;  // FIPS of the census tract to seed infected people
  unsigned int nSeedInfectedCountyFIPS; // FIPS of the county to seed infected people
  unsigned int nSeedInfectedStateFIPS;  // FIPS of the state to seed infected people
  int nSeedInfectedNumber;     // number of infected people to seed
  bool bSeedDaily;             // seed infected people every day or just once?
  int  nSeedAirports;          // seed n/10000 infected people at major airports
  bool bTravel;                // is short-term travel allowed?
  double fPreexistingImmunityFraction[TAG];  // fraction with pre-existing (sterilizing) immunity by age group
  double fPreexistingImmunityLevel; // protection against infection for those with pre-existing immunity (like VES)
  double fBaselineVESByAge[TAG];    // default VES by age group

  int nSchoolOpeningDays[56];  // day of the simulation on which the state's schools open. Values of <=0 indicate that the schools are open from the first day. Array indices correspond to the FIPS code of the state-1, starting with Alabama ("01").

  // logging parameters
  ofstream *logfile;	// the log file
  ofstream *sumfile;	// the summary file
  ofstream *individualsfile;	// the individuals file

  int nLogFileInterval;// output to log file every n days (set to 0 for no log)
  vector < unsigned int > nNumVaccineDosesUsedHistory[NUMVACCINES];// history of number of vaccine doses used (across all nodes)
  vector < unsigned int > nNumAntiviralsUsedHistory;// history of number of antiviral courses distributed (across all nodes)
  vector < unsigned int > nNumSymptomaticHistory;// history of number of symptomatics (across all nodes)
  vector < unsigned int > nNumCumulativeSymptomaticHistory;// history of number of cumulative symptomatics (across all nodes)

  // intervention parameters
  int bTrigger;                     // has the trigger for response been reached
  int nTriggerTime;                 // time when reactive strategies are deployed everywhere
  int nTriggerDelay;                // days between trigger and response
  double fResponseThreshhold;       // fraction of ever infecteds before reaction
  int nAscertainmentDelay;          // days between symptomatic and ascertainment
  double fAdultEssentialFraction;   // fraction of working-age adults who are prioritized to receive vaccine
  double fPregnantFraction[TAG];    // fraction of working-age adults who are pregnant
  double fHighRiskFraction[TAG];    // fraction of individuals who are at high risk for mortality/morbidity

  // pharmaceutical intervention parameters
  enum prevaccType ePrevaccinationStrategy; // actions to take before simulation starts
  enum vaccType eVaccinationStrategy;       // actions to take in response to epidemic
  unsigned char nVaccinePriorities[PRIORITY_LAST]; // vaccination priorities for all categories
  unsigned char nHighestPriority; // highest value in nVaccinePriority
  unsigned char nVaccinePriorities2[PRIORITY_LAST]; // vaccination priorities for all categories for "second" policy
  unsigned char nHighestPriority2; // highest value in nVaccinePriority2
  int nPriorityChangeTime;           // when to switch to "second" policy

  unsigned int nNumVaccineDosesUsed[NUMVACCINES];  // number of vaccine doses used on this node
  double fVaccinationFraction;          // vaccinate this proportion of the population
  double fContactAscertainment;         // fraction of ascertainment of non-close mixing groups
  double fSymptomaticAscertainment;    // fraction of ascertainment of symptomatic people
  unsigned int nNumAntiviralsUsed;     // number of antiviral doses used on this node
  antiviralPolicy eAVPolicy;           // use targeted antiviral prophylaxis?
  unsigned int nNumTAPDone;  // how many cases were TAPped so far.
  unsigned int nVaccineInitialSupply[NUMVACCINES];  // total doses of vaccines initially available
  unsigned int nVaccineSupply[NUMVACCINES];  // total doses of vaccines available
  unsigned int nAVTotalLimit;       // total courses of AV available
  unsigned int nVaccineDailyLimit;  // doses of vaccines that can be used per day
  unsigned int nAVDailyLimit;       // courses of AV that can be used per day
  unsigned int nNumWantVaccine;     // total number of people who want to be vaccinated
  unsigned int nNumWantAV;          // total number of people who want antivirals

  vaccinedatastruct VaccineData[NUMVACCINES];  // properties of vaccines
  int nNumVaccineTypes;                        // number of defined vaccine types
  double fVaccineEfficacyByAge[TAG];           // relative efficacy of vaccine by age
  bool bVaccineBoostByAge[TAG+3]; // requires boost?
  unsigned int vaccineproductionschedule[NUMVACCINES][MAXRUNLENGTH];
  double vload[VLOADNSUB][VLOADNDAY]; // re-scaled version of basevload from params.h

  // antiviral efficacy
  double AVEs;   // less likely to get infected
  double AVEp;   // less likely to become sick given infection
  double AVEi;   // less likely to infect others

  // non-pharmaceutical intervention parameters
  int nSchoolClosurePolicy;         // 1 for county-wide, 2 for ascertained
  int nSchoolClosureDays;           // number of days to close schools (0 for no school closure)
  double fIsolationCompliance;      // probability of voluntary home isolation compliance (set to 0 for no isolation)?
  double fQuarantineCompliance;     // probability of individual compliance (set to 0 for no quarantine)
  double fLiberalLeaveCompliance;   // probability of individual compliance (set to 0 for no liberal leave)

#ifdef PARALLEL
  int rank;                         // region(processor) number
  int size;                         // total regions(processors)
  int *tractToCPU;                  // tract id->CPU mapping
  int *numTractsPerNode;            // number of tracts on each node
  unsigned int *numPeoplePerNode;            // number of people on each node
  vector <unsigned int> emigrantworkers;     // IDs of residents working on other nodes
  vector <unsigned int> emigrantupdates;     // emigrants who were recently infected
  vector <Person> immigrantupdates; // immigrants who were recently infected
  unsigned int nNumPeopleTotal;     // number of people across all nodes
  unsigned int nNumVaccineDosesUsedTotal[NUMVACCINES];// number of vaccine doses used across all nodes
  unsigned int nNumAntiviralsUsedTotal;// number of antiviral doses used across all nodes
  MPI_Datatype PersonStubType;
  MPI_Datatype PersonStatusStubType;
#endif
};
