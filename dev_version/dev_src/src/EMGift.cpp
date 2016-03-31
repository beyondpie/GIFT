// Implementation of class EM in namespace gift.

// Libraries
#include<cmath> // for log, exp, and pow.
#include<ctime> // for c time style.
#include<chrono> // for record timing
#include<boost/thread/thread.hpp>
#include<boost/bind.hpp>
#include "gift.h"

namespace gift {

  template <typename func>
  int functionThread(func useFun,int thread, EM * point) {
    boost::thread * y;
    boost::thread_group * x = new boost::thread_group;
    for(int i=0;i<thread;++i){
      y = new boost::thread(useFun,point,i);
      x->add_thread(y);
    } // end of loop i
    x->join_all();
    delete x;
    return 0;
  } // end of function.

  double EM::iterdrugSub2ProteinSub(int drugIndex,int proteinIndex){
    double tmp = 0;
    for(auto const & m : (*drug2sub)[drugIndex]){
      for(auto const & n : (*protein2sub)[proteinIndex]){
        tmp += log(1 - (*drugSub2proteinSub)[m][n]);
      } // end of loop n
    } // end of loop m
    return exp(tmp);
  } // end of function.

  void EM::EStepThread(int threadNth){
    for(int i=threadNth;i<drugNum;i+=thread){
      for(int j=0;j<proteinNum;++j){
        double tmp = iterdrugSub2ProteinSub(i,j);
        (*observedDrug2Protein)[i][j] = (1-fn)*(1-tmp) + fp*tmp;
      } // end of for loop j
    } // end of for loop i
    //return 0;
  } // end of function

  int EM::EStep() {
    return functionThread(&EM::EStepThread,thread,this);
  } // end of function

  void EM::MStepThread(int threadNth){
    for(int i=threadNth;i<subNum;i+=thread){
      for(int j=0;j<domainNum;++j){
        double tmp = 0;
        for(auto const &m : (*sub2drug)[i]){
          for(auto const &n : (*sub2protein)[j]){
            double observed = (*observedDrug2Protein)[m][n];
            tmp += observed>0 ? (1-fn)/observed : fn/(1-observed);
          } // end of loop n
        } // end of loop m
        int tmpNum = (*sub2drug)[i].size() + (*sub2protein)[j].size();
        tmp = log((*drugSub2proteinSub)[i][j]) + log(tmp/tmpNum);
        (*drugSub2proteinSub)[i][j] = exp(tmp);
      } // end of loop j
    } // end of for loop i
  } // end of function

  int EM::MStep() {
    return functionThread(&EM::MStepThread,thread,this);
  } // end of function

  double EM::recLoglikely() {
    double loglikely = 0;
    for(int i=0;i<drugNum;++i) {
      for(int j=0;j<proteinNum;++j) {
        double tmp = iterdrugSub2ProteinSub(i,j);
        std::vector<int>::iterator it = std::find((*drug2protein)[i].begin(),
                                             (*drug2protein)[i].end(),j);
        loglikely += it==(*drug2protein)[i].end() ?
          log(1-(1-fn)*(1-tmp)-fp*tmp) : log((1-fn)*(1-tmp) + fp*tmp);
      } // end of loop j
    } // end of loop i
    return loglikely;
  } // end of function.

  int EM::trainEM(){
    // lack of loglikely record
    for(int i=0;i<iterNum;++i){
      std::cout<<"Current iteration number is " << i << std::endl;
      std::chrono::steady_clock::time_point tBegin =
        std::chrono::steady_clock::now();
      EStep();
      std::chrono::steady_clock::time_point tEnd =
        std::chrono::steady_clock::now();
      std::cout<<"Time difference (s): "
       <<std::chrono::duration_cast<std::chrono::seconds>(tBegin-tEnd).count()
               <<std::endl;

      tBegin = std::chrono::steady_clock::now();
      MStep();
      tEnd = std::chrono::steady_clock::now();
      std::cout<<"Time difference (s): "
       <<std::chrono::duration_cast<std::chrono::seconds>(tBegin-tEnd).count()
               <<std::endl;
    } // end of loop i
    std::chrono::system_clock::time_point endTime =
      std::chrono::system_clock::now();
    std::time_t endTimeT = std::chrono::system_clock::to_time_t(endTime);
    std::cout<<"Finished computation at " << std::ctime(&endTimeT) << std::endl;
    return 0;
  } // end of function

  int EM::varEM(){
    numericMatrix quesiDev;
    std::vector<double> tmpDev;
    // Note: not check drugSub2proteinsub values larger than 0.95.
    for (int i=0;i<drugNum;++i) {
      for (int j=0;j<proteinNum;++j){
        tmpDev.push_back(iterdrugSub2ProteinSub(i,j));
      } // end of loop j
      quesiDev.push_back(tmpDev);
      tmpDev.empty();
    } // end of loop i

    for(int i=0;i<subNum;++i){
      for(int j=0;j<domainNum;++j){
        double tmp_t = 0;
        //std::vector<double> tmp_t_array;
        double tmp_s = 0;
        //std::vector<double> tmp_s_array;
        double tmpLikely = (*observedDrug2Protein)[i][j];
        double tmp_sum = 0;
        for(auto const & m : (*sub2drug)[i]){
          for(auto const & n : (*sub2protein)[j]){
            tmp_t = std::find((*drug2protein)[m].begin(),(*drug2protein)[m].end(),n)
              == (*drug2protein)[m].end() ? 1/pow(1-tmpLikely,2.0) : 1/pow(tmpLikely,2.0);
            // tmp_t_array.push_back(tmp_t);
            tmp_s = (1-fn-fp) * quesiDev[m][n] / (1 - (*drugSub2proteinSub)[i][j]);
            //tmp_s_array.push_back(tmp_s);
            tmp_sum += pow(tmp_s,2.0)*tmp_t;
          } // end of loop n
        } // end of loop m
        // need initialize var matrix.
        (*vardrugSub2proteinSub)[i][j] = tmp_sum;
      } // end of loop j
    } // end of loop i
    return 0;
  } // end of function

  int EM::predictEMByDrug(IntList & drugList,
                          numericMatrix & drug2ProteinPredict){
    // Identification of existence of drugSub2proteinsub in the main function.
    //if (!drugSub2proteinSub) {
    //  std::cerr<<"ERROR: The DrugSub2ProteinSub matrix is null." <<std::endl;
    //}
    for(auto const & drug : drugList){
      for(int j=0;j<proteinNum;++j){
        drug2ProteinPredict[drug][j] = iterdrugSub2ProteinSub(drug,j);
      } // end of loop j
    } // end of loop drug
    return 0;
  } // end of function

  int EM::predictEMByProtein(IntList & proteinList,
                             numericMatrix & drug2ProteinPredict) {
    // Note: keep the matrix drug2protein.
    for(auto const & protein : proteinList){
      for(int j=0;j<drugNum;++j){
        drug2ProteinPredict[j][protein] = iterdrugSub2ProteinSub(j,protein);
      } // end of loop j
    } // end of loop protein
    return 0;
  } // end of function√

  int EM::predictEMByBoth(IntList & drugList, IntList & proteinList,
                          numericMatrix & drug2ProteinPredict) {
    for(auto const & drug : drugList){
      for(auto const & protein : proteinList){
        drug2ProteinPredict[drug][protein]=iterdrugSub2ProteinSub(drug,protein);
      } // end of loop protein
    } // end of loop drug
    return 0;
  } // end of function
} // end of namespace gift.
