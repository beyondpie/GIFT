// Implementation of class EM in namespace gift.

// Libraries
#include<cmath> // for log and exp.
#include<boost/thread/thread.hpp>
#include<boost/bind.hpp>
#include "gift.h"

namespace gift {

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
  }

  int EM::EStep() {
    boost::thread * y;
    boost::thread_group * x = new boost::thread_group;
    for (int i=0;i<thread;++i){
      y = new boost::thread(&EM::EStepThread,this,i);
      x->add_thread(y);
    }
    x->join_all();
    delete x;
    return 0;
  }

} // end of namespace gift.
