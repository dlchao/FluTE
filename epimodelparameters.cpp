/* class EpiModelParameters
 */

#include <cstdlib>
#include <cstring>
#include <climits>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "epimodel.h"
#include "epimodelparameters.h"

using namespace std;

#ifdef PARALLEL
EpiModelParameters::EpiModelParameters(int r, int s, const char *configname) {
  rank = r;
  size = s;
#else
EpiModelParameters::EpiModelParameters(const char *configname) {
#endif
  const double vacceff[VACCEFFLENGTH+1] = {0,0.001,0.005,0.015,0.033,0.061,0.1,0.152,0.218,0.301,0.401,0.519,0.658,0.818,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  // set up default values for parameters
  nRunLength = MAXRUNLENGTH;
  nLogFileInterval = 1;
  bIndividualsFile = false;
  ePrevaccinationStrategy=NOPREVACC;
  eVaccinationStrategy=NOVACC;
  for (int i=0; i<PRIORITY_LAST; i++) {
    nVaccinePriorities[i]=1;
    nVaccinePriorities2[i]=1;
  }
  nPriorityChangeTime = -1;
  fVaccinationFraction=0.7;
  for (int i=0; i<NUMVACCINES; i++) {
    nVaccineInitialSupply[i] = 0;
    memset(vaccinedata+i, 0, sizeof(vaccinedatastruct));
    memset(vaccineproductionschedule[i], 0, MAXRUNLENGTH*sizeof(unsigned int));
    for (int j=0; j<=VACCEFFLENGTH; j++)
      vaccinedata[i].vacceff[j] = vacceff[j];
  }
  fPreexistingImmunityLevel = 1.0;
  for (int i=0; i<TAG; i++) {
    fPreexistingImmunityFraction[i] = 0.0;
    fBaselineVESByAge[i] = 0.0;
    fPregnantFraction[i] = 0.0;
    fHighRiskFraction[i] = 0.0;
    fVaccineEfficacyByAge[i] = 1.0;
  }
  for (int i=0; i<TAG+3; i++) 
    bVaccineBoostByAge[i] = false;
  fResponseThreshhold=0.0;
  nAVTotalLimit = nVaccineDailyLimit = nAVDailyLimit = UINT_MAX;

  AVEs = 0.3;  // less likely to get infected
  AVEp = 0.6;  // less likely to become sick given infection
  AVEi = 0.62; // less likely to infect others

  nTriggerDelay=1;
  bTrigger=false;
  eAntiviralPolicy=NOAV;
  fContactAscertainment = 0.8;
  nSchoolClosureDays=0;
  schoolClosurePolicy=0;
  fIsolationCompliance=0.0;
  fQuarantineCompliance=0.0;
  fLiberalLeaveCompliance=0.0;
  nAscertainmentDelay = 1;
  fSymptomaticAscertainment=0.8;
  fAdultEssentialFraction = 0.0;
  nSeedInfectedTractFIPS = nSeedInfectedCountyFIPS = nSeedInfectedStateFIPS = nSeedInfectedNumber = 0;
  bSeedDaily = false;
  nSeedAirports = 0;
  bTravel = false;
  memset(nSchoolOpeningDays, 0, 56*sizeof(int));

  int flag = NoError;
  flag = (readConfigFile(configname)?flag:MiscError);
  flag = (createOutputFiles()?flag:FileError);
#ifdef PARALLEL
  MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if(flag) {
#ifdef PARALLEL
    MPI_Finalize();
#endif
    exit(flag);
  }
}

/*
 * Set up output files based on the run number read from the 'run-number' file.
 */
bool EpiModelParameters::createOutputFiles(void) {
  nr=0;
  fstream nrfile;
#ifdef PARALLEL
  if(!rank) {
#endif
    nrfile.open("run-number", ios_base::in|ios_base::out);
    if(!nrfile) {
      cerr << "WARNING: 'run-number' missing in the simulator directory. Assuming run-number is '0'." << endl;
      //      cerr << "Create a file called 'run-number' with the number 0 on the first line, then try again." << endl;
      //      return false;
      nr = 0;
    } else {
      nrfile >> nr;
      nrfile.seekg(0);
      nrfile << (nr+1) << endl;
      nrfile.close();
    }
    if (nLogFileInterval>0) {
#ifdef PARALLEL
      if (!rank) {
#else
	{
#endif
	  if (szLogFileName.length()==0) {
	    ostringstream oss;
	    oss << "Log" << nr;
	    szLogFileName = oss.str();
	  }
	EpiModelParameters::logFile.open(szLogFileName.c_str());
	if(EpiModelParameters::logFile.fail()) {
	  cerr << "ERROR: Log file '" << szLogFileName << "' cannot be open for writing." << endl;
	  return false;
	}
	if (szTractFileName.length()==0) {
	  ostringstream oss;
	  oss << "Tracts" << nr;
	  szTractFileName = oss.str();
	}
	
	EpiModelParameters::tractFile.open(szTractFileName.c_str());
	if(EpiModelParameters::tractFile.fail()) {
	  cerr << "ERROR: Tract file '" << szTractFileName << "' cannot be open for reading." << endl;
	  return false;
	}
      }
    }
#ifdef PARALLEL
  }
    MPI_Bcast(&nr, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef PARALLEL
  if (!rank) {
#else
  {
#endif
    if (szSummaryFileName.length()==0) {
      ostringstream oss;
      oss << "Summary" << nr;
      szSummaryFileName = oss.str();
    }
    summaryFile.open(szSummaryFileName.c_str());
    if(summaryFile.fail()) {
      cerr << "ERROR: Summary file '" << szSummaryFileName << "' cannot be open for writing." << endl;
      return false;
    }
  }

#ifdef PARALLEL
  if (!rank && bIndividualsFile) {
#else
  if (bIndividualsFile) {
#endif
    if (szIndividualFileName.length()==0) {
      ostringstream oss;
      oss << "Individuals" << nr;
      szIndividualFileName = oss.str();
    }
    EpiModelParameters::individualsFile.open(szIndividualFileName.c_str());
    if(EpiModelParameters::individualsFile.fail()) {
      cerr << "ERROR: Individuals file '" << szIndividualFileName << "' cannot be open for reading." << endl;
      return false;
    }
  }

  return true;
}

// read_config_bool - parses a boolean from a stream
bool EpiModelParameters::read_config_bool(bool &result, istream &iss, const char *s) {
  if (!(iss>>result)) {
#ifdef PARALLEL
    if (!rank)
#endif
      cerr << "WARNING: Invalid value for " << s << endl;
    return false;
  }
  return true;
}

// read_config_int - parses an int from a stream
bool EpiModelParameters::read_config_int(int &result, istream &iss, const char *s) {
  if (!(iss>>result)) {
#ifdef PARALLEL
    if (!rank)
#endif
      cerr << "WARNING: Invalid value for " << s << endl;
    return false;
  }
  return true;
}

// read_config_unsigned - parses an unsigned int from a stream
bool EpiModelParameters::read_config_unsigned(unsigned int &result, istream &iss, const char *s) {
  if (!(iss>>result)) {
#ifdef PARALLEL
    if (!rank)
#endif
      cerr << "WARNING: Invalid value for " << s << endl;
    return false;
  }
  return true;
}

// read_config_double - parses a double from a stream
bool EpiModelParameters::read_config_double(double &result, istream &iss, const char *s, const double min, const double max) {
  double input;
  if (!(iss>>input)) {
#ifdef PARALLEL
    if (!rank)
#endif
      cerr << "WARNING: Invalid value for " << s << endl;
    return false;
  } else if (input<min || input>max) {
#ifdef PARALLEL
    if (!rank)
#endif
      cerr << "WARNING: " << input << " is an invalid value for " << s << endl;
    return false;
  }
  result=input;
  return true;
}

/*
 * Read configuration from file `configname' to set up 
 * parameters
 */
bool EpiModelParameters::readConfigFile(const char *configname) {
  seeddisp = 0;	// should be read from config
  szBaseName = "one";
  beta     = 0.1;
  R0       = -1;
  for (int i=0; i<MAXRUNLENGTH; i++)
    seasonality[i] = 1.0;
  bool bVaccineDosesSet[NUMVACCINES];

  for (int i=0; i<NUMVACCINES; i++)
    bVaccineDosesSet[i] = false;

  // read config file
  ostringstream oss;
  if (configname)
    oss.str(configname);
  else
    oss.str("config");
  ifstream iss(oss.str().c_str());
  if (!iss) {
#ifdef PARALLEL
    if (!rank)
#endif
      cerr << "ERROR: " << configname << " not found." << endl;
#ifdef PARALLEL
    MPI_Finalize();
#endif
    exit(-1);
  } else {
    while (iss) {
      string line;
      if (getline(iss, line)) {
	istringstream iss;
	iss.str(line);
	string param;
	iss >> param;
	if (line[0]=='#') {
	  // do nothing - this line is a comment
	} else if (param.compare("label")==0) {
	  iss>>szLabel;
	} else if (param.compare("datafile")==0) {
	  iss>>szBaseName;
	} else if (param.compare("summaryfilename")==0) {
	  iss>>szSummaryFileName;
	} else if (param.compare("logfilename")==0) {
	  iss>>szLogFileName;
	} else if (param.compare("tractfilename")==0) {
	  iss>>szTractFileName;
	} else if (param.compare("individualfilename")==0) {
	  iss>>szIndividualFileName;
	} else if (param.compare("beta")==0) {
	  read_config_double(beta, iss, "beta", 0.0, 1000000.0);
	  if (beta<0.0) {
#ifdef PARALLEL
	    if (!rank)
	      cerr << "ERROR: beta can not be less than 0" << endl;
	    MPI_Finalize();
#else
	    cerr << "ERROR: beta can not be less than 0" << endl;
#endif
	    exit(-1);
	  }
	} else if (param.compare("logfile")==0) {
	  iss>>nLogFileInterval;
	} else if (param.compare("individualfile")==0) {
	  iss>>bIndividualsFile;
	} else if (param.compare("R0")==0) {
	  if (read_config_double(R0, iss, "R0", 0.0, 1000000.0)) {
	    beta = (R0-0.1099)/5.2988;
	  }
	} else if (param.compare("seasonality")==0) {
	    for (int i=0; i<MAXRUNLENGTH; i++) {
	      if (!read_config_double(seasonality[i], iss, "seasonality", 0.0, 10000.0))
		break;
	    }
	} else if (param.compare("runlength")==0) {
	    read_config_int(nRunLength, iss, "run length in days");
	    if (nRunLength>MAXRUNLENGTH) {
	      cerr << "  " << MAXRUNLENGTH << " is the maximum runlength (" << nRunLength << " is too high)." << endl;
	      exit(-1);
	    }
	} else if (param.compare("seed")==0) {
	    read_config_int(seeddisp, iss, "random seed");
	} else if (param.compare("prestrategy")==0) {
	  string s;
	  if (iss>>s) {
	    if (s.compare("prevaccinate")==0)
	      ePrevaccinationStrategy=PREVACCINATE;
	    else if (s.compare("primeboostrandom")==0)
	      ePrevaccinationStrategy=PRIMEBOOSTRANDOM;
	    else if (s.compare("primeboostsame")==0)
	      ePrevaccinationStrategy=PRIMEBOOSTSAME;
	    else if (s.compare("none")==0)
	      ePrevaccinationStrategy=NOPREVACC;
	    else {
	      cerr << "  " << s << " is not a valid pre-vaccination strategy." << endl;
	    }
	  }
	} else if (param.compare("reactivestrategy")==0) {
	  string s;
	  if (iss>>s) {
	    if (s.compare("tract")==0)
	      eVaccinationStrategy=RINGVACTRACT;
	    else if (s.compare("county")==0) {
	      eVaccinationStrategy=RINGVACCOUNTY;
	      cerr << "ERROR: County vaccination not implemented" << endl;
	    } else if (s.compare("mass")==0)
	      eVaccinationStrategy=MASSVAC;
	    else if (s.compare("none")==0)
	      eVaccinationStrategy=NOVACC;
	    else {
	      cerr << "  " << s << " is not a valid reactive vaccination strategy." << endl;
	    }
	  }
	} else if (param.compare("vaccinationfraction")==0) {
	  if (read_config_double(fVaccinationFraction, iss, "vaccination fraction", 0.0, 1.0))
	    uVaccinationFraction=(unsigned int)(UINT_MAX*fVaccinationFraction);
	} else if (param.compare("vaccinepriorities")==0) {
	  for (int i=0; i<PRIORITY_LAST; i++) {
	    int j;
	    if (!read_config_int(j, iss, "vaccine priority")) {
	      exit(-1);
	    }
	    nVaccinePriorities[i] = (unsigned char)j;
	  }
	} else if (param.compare("vaccinepriorities2")==0) {
	  for (int i=0; i<PRIORITY_LAST; i++) {
	    int j;
	    if (!read_config_int(j, iss, "second vaccine priority")) {
	      exit(-1);
	    }
	    nVaccinePriorities2[i] = (unsigned char)j;
	  }
	} else if (param.compare("prioritychangetime")==0) {
	  read_config_int(nPriorityChangeTime, iss, "priority change time");
	} else if (param.compare("antiviraldoses")==0) {
	    read_config_unsigned(nAVTotalLimit, iss, "total antiviral doses");
	} else if (param.compare("vaccinedoses")==0) {
	    int nID;
	    read_config_int(nID, iss, "vaccine ID");
	    if (nID>=NUMVACCINES) {
	      #ifdef PARALLEL
	      if (rank==0)
	      #endif
	      cerr << "ERROR: " << nID << " is not a valid vaccine ID" << endl;
	      exit(-1);
	    }
	    read_config_unsigned(nVaccineInitialSupply[nID], iss, "total initial vaccine doses");
	    bVaccineDosesSet[nID] = true;
	} else if (param.compare("antiviraldosesdaily")==0) {
	    read_config_unsigned(nAVDailyLimit, iss, "number of antiviral courses that can be administered per day");
	} else if (param.compare("vaccinedosesdaily")==0) {
	    read_config_unsigned(nVaccineDailyLimit, iss, "number of vaccine doses that can be administered per day");
	} else if (param.compare("vaccinedata")==0) {
	    int nID;
	    read_config_int(nID, iss, "vaccine ID");
	    if (nID>=NUMVACCINES) {
	      #ifdef PARALLEL
	      if (rank==0)
	      #endif
	      cerr << "ERROR: " << nID << " is not a valid vaccine ID" << endl;
	      exit(-1);
	    }
	    if (vaccinedata[nID].nNumDoses<=0)
	      vaccinedata[nID].nNumDoses = 1;
	    if (!bVaccineDosesSet[nID])
	      nVaccineInitialSupply[nID]=UINT_MAX;
	    read_config_double(vaccinedata[nID].VEs, iss, "VEs", 0.0, 1.0);
	    read_config_double(vaccinedata[nID].VEi, iss, "VEi", 0.0, 1.0);
	    read_config_double(vaccinedata[nID].VEp, iss, "VEp", 0.0, 1.0);
	    read_config_double(vaccinedata[nID].fInfants, iss, "bad for infants", 0.0, 1.0);
	    for (int i=0; i<TAG; i++)
	      read_config_double(vaccinedata[nID].fAge[i], iss, "bad for age group i", 0.0, 1.0);
	    if (!read_config_bool(vaccinedata[nID].bPregnant, iss, "bad for pregnant women")) {
	      
	      #ifdef PARALLEL
	      if (rank==0)
	      #endif
	      cerr << "ERROR: error in vaccinedata values for vaccine " << nID << endl;	      
	      exit(-1);
	    }
	} else if (param.compare("vaccineproduction")==0) {
	    int nID;
	    read_config_int(nID, iss, "vaccine ID");
	    if (nID>=NUMVACCINES) {
	      #ifdef PARALLEL
	      if (rank==0)
	      #endif
	      cerr << "ERROR: " << nID << " is not a valid vaccine ID" << endl;
	      exit(-1);
	    }
	    for (int i=0; i<MAXRUNLENGTH; i++) {
	      if (!read_config_unsigned(vaccineproductionschedule[nID][i], iss, "daily vaccine production"))
		break;
	    }
	} else if (param.compare("AVEs")==0) {
	    read_config_double(AVEs, iss, "AVEs", 0.0, 1.0);
	} else if (param.compare("AVEi")==0) {
	    read_config_double(AVEi, iss, "AVEi", 0.0, 1.0);
	} else if (param.compare("AVEp")==0) {
	    read_config_double(AVEp, iss, "AVEp", 0.0, 1.0);
	} else if (param.compare("responsethreshhold")==0) {
	    read_config_double(fResponseThreshhold, iss, "response threshhold", 0.0, 1.0);
	} else if (param.compare("responsedelay")==0) {
	    read_config_int(nTriggerDelay, iss, "response delay");
	} else if (param.compare("responseday")==0) {
	    read_config_int(nTriggerDay, iss, "response day");
	} else if (param.compare("ascertainmentdelay")==0) {
	    read_config_int(nAscertainmentDelay, iss, "ascertainment delay");
	} else if (param.compare("ascertainmentfraction")==0) {
	    read_config_double(fSymptomaticAscertainment, iss, "ascertainment fraction", 0.0, 1.0);
	} else if (param.compare("essentialfraction")==0) {
	    read_config_double(fAdultEssentialFraction, iss, "essential workforce fraction", 0.0, 1.0);
	} else if (param.compare("pregnantfraction")==0) {
	  read_config_double(fPregnantFraction[0], iss, "preschool pregnant fraction", 0.0, 1.0);
	  read_config_double(fPregnantFraction[1], iss, "school-age pregnant fraction", 0.0, 1.0);
	  read_config_double(fPregnantFraction[2], iss, "young adult pregnant fraction", 0.0, 1.0);
	  read_config_double(fPregnantFraction[3], iss, "older adult pregnant fraction", 0.0, 1.0);
	  read_config_double(fPregnantFraction[4], iss, "elderly pregnant fraction", 0.0, 1.0);
	} else if (param.compare("highriskfraction")==0) {
	    read_config_double(fHighRiskFraction[0], iss, "high risk fraction", 0.0, 1.0);
	    read_config_double(fHighRiskFraction[1], iss, "high risk fraction", 0.0, 1.0);
	    read_config_double(fHighRiskFraction[2], iss, "high risk fraction", 0.0, 1.0);
	    read_config_double(fHighRiskFraction[3], iss, "high risk fraction", 0.0, 1.0);
	    read_config_double(fHighRiskFraction[4], iss, "high risk fraction", 0.0, 1.0);
	} else if (param.compare("seedtract")==0) {
	  if (iss>>nSeedInfectedStateFIPS && 
	      iss>>nSeedInfectedCountyFIPS &&
	      iss>>nSeedInfectedTractFIPS &&
	      iss>>nSeedInfectedNumber) {
	  } else {
	    cerr << "seedtract failed" << endl;
	  }
	} else if (param.compare("seedinfected")==0) {
	    if (!(iss>>nSeedInfectedNumber)) {
	      cerr << "seedinfected failed" << endl;
	    }
	} else if (param.compare("seedinfecteddaily")==0) {
	    read_config_bool(bSeedDaily, iss, "daily seeding");
	} else if (param.compare("seedairports")==0) {
	    read_config_int(nSeedAirports, iss, "airport seeding");
	} else if (param.compare("travel")==0) {
	    read_config_bool(bTravel, iss, "travel");
	} else if (param.compare("antiviralpolicy")==0) {
	  string s;
	  if (iss>>s) {
	    if (s.compare("none")==0)
	      eAntiviralPolicy=NOAV;
	    else if (s.compare("treatmentonly")==0)
	      eAntiviralPolicy=TREATMENTONLY;
	    else if (s.compare("HHTAP")==0)
	      eAntiviralPolicy=HHTAP;
	    else if (s.compare("HHTAP100")==0)
	      eAntiviralPolicy=HHTAP100;
	    else if (s.compare("FULLTAP")==0)
	      eAntiviralPolicy=FULLTAP;
	    else {
	      cerr << "  " << s << " is not a valid antiviral strategy." << endl;
	    }
	  }
	} else if (param.compare("schoolopeningdays")==0) {
	  for (int i=0; i<56; i++)
	    if (!read_config_int(nSchoolOpeningDays[i], iss, "school opening days"))
	      exit(-1);
	} else if (param.compare("schoolclosuredays")==0) {
	    read_config_int(nSchoolClosureDays, iss, "school closure days");
	} else if (param.compare("schoolclosurepolicy")==0) {
	  string s;
	    if (iss>>s) {
	      if (s.compare("none")==0)
		schoolClosurePolicy=0;
	      else if (s.compare("all")==0)
		schoolClosurePolicy=1;
	      else if (s.compare("bytractandage")==0)
		schoolClosurePolicy=2;
	      else {
		cerr << "WARNING: " << s << " is not a valid school closure strategy." << endl;
	      }
	    }
	} else if (param.compare("isolation")==0) {
	    read_config_double(fIsolationCompliance, iss, "isolation compliance probability", 0.0, 1.0);
	} else if (param.compare("quarantine")==0) {
	    read_config_double(fQuarantineCompliance, iss, "quarantine compliance probability", 0.0, 1.0);
	} else if (param.compare("liberalleave")==0) {
	    read_config_double(fLiberalLeaveCompliance, iss, "liberal leave compliance probability", 0.0, 1.0);
	} else if (param.compare("preexistingimmunitylevel")==0) {
	    read_config_double(fPreexistingImmunityLevel, iss, "protection for those with pre-existing immunity", 0.0, 1.0);
	} else if (param.compare("preexistingimmunitybyage")==0) {
	  read_config_double(fPreexistingImmunityFraction[0], iss, "fraction of pre-schoolers with pre-existing immunity", 0.0, 1.0);
	  read_config_double(fPreexistingImmunityFraction[1], iss, "fraction of school-age children with pre-existing immunity", 0.0, 1.0);
	  read_config_double(fPreexistingImmunityFraction[2], iss, "fraction of young adults with pre-existing immunity", 0.0, 1.0);
	  read_config_double(fPreexistingImmunityFraction[3], iss, "fraction of older adults with pre-existing immunity", 0.0, 1.0);
	  read_config_double(fPreexistingImmunityFraction[4], iss, "fraction of elderly with pre-existing immunity", 0.0, 1.0);
	} else if (param.compare("defaultVESbyage")==0) {
	  read_config_double(fBaselineVESByAge[0], iss, "default VES for pre-schoolers", 0.0, 1.0);
	  read_config_double(fBaselineVESByAge[1], iss, "default VES for school-age children", 0.0, 1.0);
	  read_config_double(fBaselineVESByAge[2], iss, "default VES for young adults", 0.0, 1.0);
	  read_config_double(fBaselineVESByAge[3], iss, "default VES for older adults", 0.0, 1.0);
	  read_config_double(fBaselineVESByAge[4], iss, "default VES for elderly", 0.0, 1.0);
	} else if (param.compare("vaccineefficacybyage")==0) {
	  read_config_double(fVaccineEfficacyByAge[0], iss, "relative vaccine efficacy for pre-schoolers", 0.0, 1.0);
	  read_config_double(fVaccineEfficacyByAge[1], iss, "relative vaccine efficacy for school-age children", 0.0, 1.0);
	  read_config_double(fVaccineEfficacyByAge[2], iss, "relative vaccine efficacy for young adults", 0.0, 1.0);
	  read_config_double(fVaccineEfficacyByAge[3], iss, "relative vaccine efficacy for older adults", 0.0, 1.0);
	  read_config_double(fVaccineEfficacyByAge[4], iss, "relative vaccine efficacy for elderly", 0.0, 1.0);
	} else if (param.compare("vaccineboostbyage")==0) {
	  read_config_bool(bVaccineBoostByAge[0], iss, "preschoolers need boosts");
	  read_config_bool(bVaccineBoostByAge[1], iss, "school-age children need boosts");
	  read_config_bool(bVaccineBoostByAge[2], iss, "young adults need boosts");
	  read_config_bool(bVaccineBoostByAge[3], iss, "older adults need boosts");
	  read_config_bool(bVaccineBoostByAge[4], iss, "elderly need boosts");
	  read_config_bool(bVaccineBoostByAge[5], iss, "elementary schoolers need boosts");
	  read_config_bool(bVaccineBoostByAge[6], iss, "middle schoolers need boosts");
	  read_config_bool(bVaccineBoostByAge[7], iss, "high schoolers need boosts");
	} else if (param.compare("vaccinebuildup")==0) {
	    int nID;
	    read_config_int(nID, iss, "vaccine ID");
	    if (nID>=NUMVACCINES) {
	      cerr << "ERROR: " << nID << " is not a valid vaccine ID" << endl;
	      exit(-1);
	    }
	    read_config_int(vaccinedata[nID].nBoostDay, iss, "vaccine boost day");
	    if (vaccinedata[nID].nBoostDay>0)
	      vaccinedata[nID].nNumDoses=2;
	    else
	      vaccinedata[nID].nNumDoses=1;
	    for (int i=0; i<=VACCEFFLENGTH; i++) {
	      if (iss>>vaccinedata[nID].vacceff[i]) {
		if (vaccinedata[nID].vacceff[i]>1.0 || vaccinedata[nID].vacceff[i]<0.0) {
		  cerr << "ERROR: invalid vaccine efficacy for vaccine " << nID << " for day " << i << "." << endl;

		  exit(-1);
		}
	      } else {
		vaccinedata[nID].vacceff[i] = 1.0;
		cerr << "WARNING: Vaccine efficacy is not defined on day " << i << " for vaccine " << nID << "." << endl;
	      }
	    }
	} else if (!param.compare("")==0) {
#ifdef PARALLEL
	  if (!rank)
#endif
	    cerr << "WARNING: " << param << " is not a valid parameter" << endl;
	}
      }
    }
    iss.close();
  }
  return true;
}
