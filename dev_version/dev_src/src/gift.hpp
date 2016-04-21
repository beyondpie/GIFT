// The head file for the method of gift.

// GIFT is used to predicte the interactions between compounds
// substructures and protein domains from the known drug-protein
// interactions. See the details:
// http://bioinformatics.oxfordjournals.org/content/31/15/2523.abstract

// Note: in this version, we don't support
//   1. the association method.
//   2. the cross-validation method and draw the AUC curve.
//   3. check for drugs/protein without fingerprints.
//   4. check for fingerpints involved in few drugs/proteins.

// Author: Songpeng Zu
// Email: zusongpeng@gmail.com
// Date: 2016-03-06

#ifndef __GIFT_H__
#define __GIFT_H__

// C/C++ system library
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
// Boost library
#include<boost/thread/thread.hpp>
#include<boost/bind.hpp>
//#include<boost/thread/thread.hpp>
namespace gift {
  // gift information.
  const std::string author("Songpeng Zu");
  const std::string email("zusongpeng@gmail.com");
  const std::string version("gift-2.0");
  const std::string updateTime("2016-03-06");
  const int recLogLeastNum = 5;

  typedef std::vector<int> IntList;
  typedef std::vector<std::vector<int> > IntArrayList;
  typedef std::vector<std::vector<double> > numericMatrix;
  typedef std::map<std::string,int> name2IndexHash;
  typedef std::vector<std::string> nameList;

  extern IntArrayList drug2proteinList;
  extern IntArrayList drug2subList;
  extern IntArrayList sub2drugList;
  extern IntArrayList protein2domainList;
  extern IntArrayList domain2proteinList;
  extern numericMatrix drugSub2proteinSubMatrix;
  extern numericMatrix observedDrug2ProteinMatrix;
  extern numericMatrix vardrugSub2proteinSubMatrix;
  extern std::vector<double> loglikelyArray;

  extern name2IndexHash drugName2Index;
  extern name2IndexHash proteinName2Index;
  extern nameList drugNameList;
  extern nameList proteinNameList;
  extern nameList drugSubNameList;
  extern nameList proteinSubNameList;

  extern nameList predictDrugNameList;
  extern nameList predictProteinNameList;

  extern nameList predictDrugNameList_WithSubs;
  extern nameList predictProteinNameList_WithSubs;
  extern IntArrayList predictDrug2SubList;
  extern IntArrayList predictProtein2SubList;

  class rowCol;
  class parameters;
  class EM;

  // gift global functions.
  int Matrix2Fingerpints(const std::string, IntArrayList&,
                         std::string delims="\t,");
  int Matrix2FingerprintsByColumn(const std::string, IntArrayList&, int rowNum,
                                  std::string delims="\t,");
  int writeMatrix(const std::string, numericMatrix&, std::string delims="\t,");
  int readMatrix(const std::string, numericMatrix&, std::string delims="\t,");

  int readNameListFromFile(const std::string, nameList&);
  int readNameMatrixFromFile(const std::string, nameList&, IntArrayList&,
                             std::string delims="\t,");

  int readName2IndexHash(const nameList, name2IndexHash&);
  int getIndexFromHash(const name2IndexHash&, const nameList, IntList&,
                       nameList&);

  int rowColFile(const std::string, rowCol&, std::string delims="\t,");
  int helpGift();
  int outRecord(parameters&, EM&);

  // put template or inline function in one file.
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

  // classes definition.
  class rowCol {
  public:
    rowCol(int row, int col): rowNum(row), colNum(col) {}
    rowCol(): rowNum(1), colNum(1) {}
    int rowNum;
    int colNum;
  };
} // end of namepsace gift

#endif // end of GIFT_H
