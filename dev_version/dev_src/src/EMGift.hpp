// head file for EMGift.
#ifndef __EMGIFT_H__
#define __EMGIFT_H__
// Libraries
#include<cmath> // for log, exp, and pow.
#include<ctime> // for c time style.
#include<chrono> // for record time
#include<boost/algorithm/string.hpp>
#include<boost/algorithm/string/join.hpp>
#include<boost/range/adaptor/transformed.hpp>
#include "gift.hpp"
#include "ParameterGift.hpp"

namespace gift{
  class EM {
    // In fact, usually only one EM objact is allowed.
  public:
    // Inition with parameters
    EM(parameters& param)
      : loglikelyRecord(param.loglikelyRecord)
      , fn(param.fn)
      , fp(param.fp)
      , thread(param.thread)
      , iterNum(param.iterNum)
      , drugNum(param.drugNum)
      , proteinNum(param.proteinNum)
      , task(param.task)
      , drug2sub(&drug2subList)
      , sub2drug(&sub2drugList)
      , protein2sub(&protein2domainList)
      , sub2protein(&domain2proteinList)
      , drug2protein(&drug2proteinList)
      , drugSub2proteinSub(&drugSub2proteinSubMatrix)
      , observedDrug2Protein(&observedDrug2ProteinMatrix)
      , vardrugSub2proteinSub(&vardrugSub2proteinSubMatrix)
      , loglikely(&loglikelyArray)
      , predictDrugsFileName(param.predictDrugsFileName)
      , predictProteinsFileName(param.predictProteinsFileName)
      , predictDrugsFileName_WithSubs(param.predictDrugsFileName_WithSubs)
      , predictProteinsFileName_WithSubs(param.predictProteinsFileName_WithSubs)
      , outputDelims(param.outputDelims)
      , outRecordFileName(param.outRecordFileName)
      , outPredictCPIsFileName(param.outPredictCPIsFileName)
      , outDrugSub2ProteinSubFileName(param.outDrugSub2ProteinSubFileName)
      , outVarDrugSub2proteinSubFileName(param.outVarDrugSub2proteinSubFileName)
    { } // end of constuctor.

    ~EM() { } // end of default destruction.
    // Set the pointers to several matrix.
    inline int setPointerDrug2Sub(IntArrayList & d2s) {
      drug2sub = &d2s;
      return 0;
    } // end of func
    inline int setPointerProtein2Sub(IntArrayList & p2s) {
      protein2sub = &p2s;
      return 0;
    } // end of func
    inline int setPointerDrug2Protein(IntArrayList & d2p) {
      drug2protein = &d2p;
      return 0;
    } // end of func
    inline int setPointerDrugSub2ProteinSub(numericMatrix & ds2ps) {
      drugSub2proteinSub = &ds2ps;
      return 0;
    } // end of func

    // two init EM, i.e., the drugSub2proteinSub matrix.
    // when the task is train, we use int initEM() to init the matrix
    // when the task is predict, we directoly read the trained matrix.
    //    int initEM();
    // Core of EM.
    double iterdrugSub2ProteinSub(int drugIndex, int proteinIndex);
    void EStepThread(int threadNth);
    int EStep();
    void MStepThread(int threadNth);
    int MStep();
    double recLoglikely();
    inline int setLoglikely(double logscore) {
      (*loglikely).push_back(logscore);
      return 0;
    } // end of function
    int trainEM();
    // predictEM.
    // int predictEMByDrug(IntList &, numericMatrix &);
    // int predictEMByProtein(IntList &, numericMatrix &);
    // int predictEMByBoth(IntList & drugs, IntList & proteins, numericMatrix &);
    int predictDrugs();
    int predictProteins();
    int predictDrugsWithSubs();
    int predictProteinsWithSubs();
    int predictDrugsWithSubsProteinsWithSubs();
    int predictDrugsProteinsWithSubs();
    int predictDrugsWithSubsProteins();
    int predictDrugsProteins();

    int predictEM();
    int varEM();
    int outTrainResult();
    int outTrainVariance();
    //int outPredict(std::string);
  private:
    bool loglikelyRecord;
    double fn;
    double fp;
    int thread;
    int iterNum;
    int drugNum;
    int subNum;
    int domainNum;
    int proteinNum;
    std::string task;
    // pointer to outside data.
    IntArrayList * drug2sub;
    IntArrayList * sub2drug;
    IntArrayList * protein2sub;
    IntArrayList * sub2protein;
    IntArrayList * drug2protein;
    numericMatrix * drugSub2proteinSub;
    numericMatrix * observedDrug2Protein;
    numericMatrix * vardrugSub2proteinSub;
    std::vector<double> * loglikely;

    // input file names for prediction.
    std::string predictDrugsFileName;
    std::string predictProteinsFileName;
    std::string predictDrugsFileName_WithSubs;
    std::string predictProteinsFileName_WithSubs;

    // output file name and format
    std::string outputDelims;
    std::string outRecordFileName;
    std::string outPredictCPIsFileName;
    std::string outDrugSub2ProteinSubFileName;
    std::string outVarDrugSub2proteinSubFileName;
  };// end of class EM

} // end of namespace gift.

#endif // endif for __EMGIFT_H__
