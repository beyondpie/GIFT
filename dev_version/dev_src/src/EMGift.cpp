// Implementation of class EM in namespace gift.

// Libraries
#include<cmath> // for log and exp.
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

  void EM::EStepThread(int threadNth){
    for(int i=threadNth;i<drugNum;i+=thread){
      for(int j=0;j<proteinNum;++j){
        double tmp = 0;
        for (auto const &m : (*drug2sub)[i]) {
          for (auto const &n : (*protein2sub)[j]) {
            tmp += log( 1 - (*drugSub2proteinSub)[m][n] );
          } // end of for loop n
        } // end of for loop m
        tmp = exp(tmp);
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
        double tmp = 0;
        for(auto const & m : (*drug2sub)[i]) {
          for(auto const & n : (*protein2sub)[j]) {
            tmp += log( 1 - (*drugSub2proteinSub)[m][n]);
          } // end of loop n
        } // end of loop m
        tmp = exp(tmp);
        std::vector<int>::iterator it = std::find((*drug2protein)[i].begin(),
                                             (*drug2protein)[i].end(),j);
        loglikely += it==(*drug2protein)[i].end() ?
          log(1-(1-fn)*(1-tmp)-fp*tmp) : log((1-fn)*(1-tmp) + fp*tmp);
      } // end of loop j
    } // end of loop i
    return loglikely;
  } // end of function.

  int EM::trainEM(){
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

} // end of namespace gift.
