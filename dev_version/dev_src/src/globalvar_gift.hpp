// global variables for gift.
#ifndef __GLOBALVAR_GIFT__
#define __GLOBALVAR_GIFT__

#include "gift.hpp"
namespace gift{
  // gift global variable.
  IntArrayList drug2proteinList;
  IntArrayList drug2subList;
  IntArrayList sub2drugList;
  IntArrayList protein2domainList;
  IntArrayList domain2proteinList;
  numericMatrix drugSub2proteinSubMatrix;
  numericMatrix observedDrug2ProteinMatrix;
  numericMatrix vardrugSub2proteinSubMatrix;
  std::vector<double> loglikelyArray;

  name2IndexHash drugName2Index;
  name2IndexHash proteinName2Index;
  nameList drugNameList;
  nameList proteinNameList;
  nameList drugSubNameList;
  nameList proteinSubNameList;

  nameList predictDrugNameList;
  nameList predictProteinNameList;

  nameList predictDrugNameList_WithSubs;
  nameList predictProteinNameList_WithSubs;
  IntArrayList predictDrug2SubList;
  IntArrayList predictProtein2SubList;

}

#endif // __GLOBALVAR_GIFT__
