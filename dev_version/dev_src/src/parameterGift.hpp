// head file for class parameters.
#ifndef __PARAMETERGIFT_H__
#define __PARAMETERGIFT_H__

#include<boost/program_options.hpp>
#include<boost/any.hpp>
#include "gift.hpp"

namespace gift {
  class parameters {
  // Further design recommmendation:
  // parameters should be the SINGLETON, a specific desing pattern.
  // And the outside in gift  data should be stored inside parameters
  //  with shared ptr.
  // Note: some people did not agree with traditional singleton, since it
  // is just like global variables. And it is not easy to keep thread safe
  // since the use of static variables.
  public:
    //parameters ();
    // Init with config file.
    parameters (const std::string) throw(std::string);

    inline int setDrugNum (int number) { drugNum = number; return 0; }
    inline int setSubNum (int number) { subNum = number; return 0; }
    inline int setProteinNum (int number) {proteinNum = number; return 0; }
    inline int setDomainNum (int number) { domainNum = number; return 0; }

    int InitDrugSub2ProteinSub();
    int InitVarDrugSub2ProteinSub();
    int InitDrugName2Index();
    int InitProteinName2Index();
    int InitDrugSubNameList();
    int InitProteinSubNameList();
    int InitPredictParameters() throw(std::string);

    // DATA MEMBERS
    // input data file name
    std::string drug2proteinFileName;
    std::string drug2subFileName;
    std::string protein2subFileName;
    std::string drugSub2proteinSubFileName;
    // input name list file name
    std::string drugNameListFile;
    std::string drugSubNameListFile;
    std::string proteinNameListFile;
    std::string proteinSubNameListFile;
    // input parameters for EM
    bool loglikelyRecord;
    double alphaEB; // Emprical Bayesian estimated parameter.
    double betaEB; // Emprical Bayesian estimated parameter.
    double fn;
    double fp;
    int thread;
    int iterNum;
    int drugNum;
    int subNum;
    int domainNum;
    int proteinNum;
    std::string inputDelims;
    std::string task;
    // input file version information.
    std::string chemfpRec;
    std::string proteinfpRec;
    std::string CPIsRec;
    // input file names for prediction.
    // If both files are given, we only predict the interactions between them.
    // It means we now support only one-time prediction.
    std::string predictDrugsFileName;
    std::string predictProteinsFileName;
    // files with subs, each line is a compound/protein, and first column is the
    // compound/protein Names.
    std::string predictDrugsFileName_WithSubs;
    std::string predictProteinsFileName_WithSubs;
    // output file name and format
    std::string outputDelims;
    std::string outRecordFileName;
    std::string outPredictCPIsFileName;
    std::string outDrugSub2ProteinSubFileName;
    std::string outVarDrugSub2proteinSubFileName;
  }; // end of class parameters
} // end of namespace gift.

#endif // end of __PARAMETERGIFT_H__
