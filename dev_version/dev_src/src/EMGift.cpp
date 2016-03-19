// Implementation of class EM in namespace gift.

// Libraries
#include "gift.h"

namespace gift {
  EM::EM(parameters& param)
    : fn(param::fn)
    , fp(param::fp)
    , thread(param::thread)
    , iterNum(param::iterNum)
    , drugNum(param::drugNum)
    , proteinNum(param::proteinNum)
    , task(param::task) {
    drug2sub = nullptr; // c++11
    protein2sub = nullptr;
    drug2protein = nullptr;
    drugSub2proteinSub = nullptr;
    loglikely = nullptr;
  } // end of class EM constructor.

} // end of namespace gift.
