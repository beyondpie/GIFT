// Implementation of class EM in namespace gift.
#include<cmath> // for log, exp, and pow.
#include<ctime> // for c time style.
#include<chrono> // for record time.
#include<boost/algorithm/string.hpp>
#include<boost/algorithm/string/join.hpp>
#include<boost/range/adaptor/transformed.hpp>

#include "gift.hpp"

namespace gift {
  double EM::iterdrugSub2ProteinSub(int drugIndex,int proteinIndex){
    double tmp = 0;
    if ((*drug2sub)[drugIndex].empty() || (*protein2sub)[proteinIndex].empty()){
      std::cerr<<"[ERROR]: some row of drug2sub or protein2sub is empty..."<<std::endl;
      return 1;
    } // end of if
    for(auto const & m : (*drug2sub)[drugIndex]){
      for(auto const & n : (*protein2sub)[proteinIndex]){
        tmp += log(1 - (*drugSub2proteinSub)[m][n]);
      } // end of loop n
    } // end of loop m
    return exp(tmp);
  } // end of function.

  void EM::EStepThread(int threadNth){
    //std::cout<<"This is thread "<<threadNth<<" for EStep..."<<std::endl;
    for(int i=threadNth;i<drugNum;i+=thread){
      for(int j=0;j<proteinNum;++j){
        double tmp = iterdrugSub2ProteinSub(i,j);
        (*observedDrug2Protein).at(i).at(j) = (1-fn)*(1-tmp) + fp*tmp;
      } // end of for loop j
    } // end of for loop i
    //return 0;
  } // end of function

  int EM::EStep(int){
    for(int i=0;i<drugNum;++i){
      for(int j=0;j<proteinNum;++j){
        double tmp = iterdrugSub2ProteinSub(i,j);
        observedDrug2ProteinMatrix[i][j] = (1-fn)*(1-tmp) + fp*tmp;
        //(*observedDrug2Protein).at(i).at(j) = (1-fn)*(1-tmp) + fp*tmp;
      } // end of for loop j
    } // end of for loop i
    return 0;
  } // end of function

  int EM::EStep() {
    // both of them works.
    return gift::functionThread(&EM::EStepThread,thread,this);
    //return functionThread(&EM::EStepThread, thread);
  } // end of function

  void EM::MStepThread(int threadNth){
    //std::cout<<"This is thread "<<threadNth<<" for MStep..."<<std::endl;
    for(int i=threadNth;i<subNum;i+=thread){
      for(int j=0;j<domainNum;++j){
        double tmp = 0;
        for(auto const &m : (*sub2drug)[i]){
          for(auto const &n : (*sub2protein)[j]){
            double observed = (*observedDrug2Protein)[m][n];
            tmp = (std::find((*drug2protein)[m].begin(),(*drug2protein)[m].end(),n)
              != (*drug2protein)[m].end() ) ? (1-fn)/observed : fn/(1-observed);
          } // end of loop n
        } // end of loop m
        int tmpNum = (*sub2drug)[i].size() + (*sub2protein)[j].size();
        tmp = log((*drugSub2proteinSub)[i][j]) + log(tmp/tmpNum);
        (*drugSub2proteinSub)[i][j] = exp(tmp);
      } // end of loop j
    } // end of for loop i
  } // end of function

  int EM::MStep(int){
    for(int i=0;i<subNum;++i){
      for(int j=0;j<domainNum;++j){
        double tmp = 0;
        for(auto  &m : (*sub2drug)[i]){
          for(auto  &n : (*sub2protein)[j]){
            double observed = (*observedDrug2Protein)[m][n];
            tmp = (std::find((*drug2protein)[m].begin(),(*drug2protein)[m].end(),n)
              != (*drug2protein)[m].end() ) ? (1-fn)/observed : fn/(1-observed);
          } // end of loop n
        } // end of loop m
        int tmpNum = (*sub2drug)[i].size() + (*sub2protein)[j].size();
        tmp = log((*drugSub2proteinSub)[i][j]) + log(tmp/tmpNum);
        drugSub2proteinSubMatrix[i][j] = exp(tmp);
        //        (*drugSub2proteinSub).at(i).at(j) = exp(tmp);
      } // end of loop j
    } // end of for loop i
    return 0;
  } // end of function

  int EM::MStep() {
    // both of them works.
    return gift::functionThread(&EM::MStepThread,thread,this);
    //return functionThread(&EM::MStepThread,thread);
  } // end of function

  double EM::recLoglikely() {
    double loglikely = 0;
    double tmp;
    std::vector<int>::iterator it;
    for(int i=0;i<drugNum;++i) {
      for(int j=0;j<proteinNum;++j) {
        tmp = iterdrugSub2ProteinSub(i,j);
        it = std::find((*drug2protein)[i].begin(),
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
      //std::cout<<"EStep Testing..."<<std::endl;
      //EStep(1); // for testEM
      std::chrono::steady_clock::time_point tEnd =
        std::chrono::steady_clock::now();
      std::cout<<"EStep Time difference (s): "
       <<std::chrono::duration_cast<std::chrono::seconds>(tBegin-tEnd).count()
               <<std::endl;

      tBegin = std::chrono::steady_clock::now();
      MStep();
      //MStep(1); // for testEM
      tEnd = std::chrono::steady_clock::now();
      std::cout<<"MStep Time difference (s): "
       <<std::chrono::duration_cast<std::chrono::seconds>(tBegin-tEnd).count()
               <<std::endl;

      if (loglikelyRecord || i<=recLogLeastNum || i >= iterNum-recLogLeastNum) {
        double tmplog = recLoglikely();
        std::cout<<"Current loglikelyhood is " << tmplog << std::endl;
        // // for test
        // std::cout<<"Current observedDrug2Protein Matrix is: "<<std::endl;
        // printMatrix(*observedDrug2Protein);
        // std::cout<<"Current drugSub2proteinSub Matrix is: "<<std::endl;
        // printMatrix(*drugSub2proteinSub);

        setLoglikely(tmplog);
      } // end of if
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
      tmpDev.clear();
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
        (*vardrugSub2proteinSub)[i][j] = tmp_sum;
      } // end of loop j
    } // end of loop i
    return 0;
  } // end of function

  int EM::predictEM(){
    // Note: we only run one situation one time:
    // - Provide both drugs and proteins Names.
    // - Provide both drugs and proteins with Subs.
    // - Provide both drugs Names and proteins with Subs.
    // - Provide both drugs with Subs and proteins Names.
    // - Provide only drugs Names.
    // - Provide only drug with Subs.
    // - Provide only protein Names.
    // - Provide only protein with Subs.
    std::cout<<"Now Run predictEM for task: predict..." << std::endl;

    if (!predictDrugNameList_WithSubs.empty() &
        !predictProteinNameList_WithSubs.empty()){
      predictDrugsWithSubsProteinsWithSubs();
    }  else if (!predictDrugNameList_WithSubs.empty() &
        !predictProteinNameList.empty()){
      predictDrugsWithSubsProteins();
    } else if (!predictProteinNameList_WithSubs.empty() &
        !predictDrugNameList.empty()){
      predictDrugsProteinsWithSubs();
    } else if (!predictDrugNameList.empty() &
        !predictProteinNameList.empty()) {
      predictDrugsProteins();
    } else if (!predictDrugNameList.empty()){
      predictDrugs();
    } else if (!predictDrugNameList_WithSubs.empty()){
      predictDrugsWithSubs();
    } else if (!predictProteinNameList.empty()){
      predictProteins();
    } else if (!predictProteinNameList_WithSubs.empty()){
      predictProteinsWithSubs();
    } else {
      std::cerr<<"No files for task: predict, and quit."<< std::endl;
      return 1;
    } // end of if else if else.
    return 0;
  } // end of function

  int EM::predictDrugs(){
    std::cout<<"Now predict given drugs against all the proteins, "
             <<"which are in our training data set." << std::endl;
    IntList predictDrugIndex;
    std::vector<double> tmpCalc;
    nameList existNameList;
    getIndexFromHash(drugName2Index, predictDrugNameList,
                     predictDrugIndex, existNameList);
    if (existNameList.size() < 1) {
      std::cerr<<"No drug Index found, and quit." <<std::endl;
      return 1;
    } // end of if
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as protein names.
    output<<"Names"<<outputDelims;
    output<<join(proteinNameList,outputDelims)<<std::endl;
    int num = 0;
    for(const auto & drug : predictDrugIndex){
      for(int j=0;j<proteinNum;++j){
        tmpCalc.push_back(iterdrugSub2ProteinSub(drug,j));
      } // end of loop for j
      output<<existNameList[num] << outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
      ++num;
    } // end of loop for drug
    return 0;
  } // end of functions.

  int EM::predictDrugsWithSubs(){
    std::cout<<"Now predict given drugs with subs against all the proteins, "
             <<"which are in our training data set." << std::endl;
    std::vector<double> tmpCalc;
    nameList existNameList;
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as protein names.
    output<<"Names"<<outputDelims;
    output<<join(proteinNameList,outputDelims)<<std::endl;
    double tmp;
    for(int i=0;i<predictDrugNameList_WithSubs.size();++i){
      for(int j=0;j<proteinNum;++j){
        tmp = 0;
        for(auto const & m : predictDrug2SubList[i]){
          for(auto const & n : (*protein2sub)[j]){
            tmp += log(1 - (*drugSub2proteinSub)[m][n]);
          } // end of loop for n
        } // end of loop for m
        tmpCalc.push_back(exp(tmp));
      } // end of loop for j
      output<<predictDrugNameList_WithSubs[i]<<outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
    } // end of loop for drug
    return 0;
  } // end of functions.

  int EM::predictProteins(){
    std::cout<<"Now predict given proteins against all the drugs, "
             <<"which are in our training data set." << std::endl;
    IntList predictProteinsIndex;
    std::vector<double> tmpCalc;
    nameList existNameList;
    getIndexFromHash(proteinName2Index, predictProteinNameList,
                     predictProteinsIndex, existNameList);
    if (existNameList.size() < 1) {
      std::cerr<<"No protein Index found, and quit." <<std::endl;
      return 1;
    } // end of if
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as drug names.
    output<<"Names"<<outputDelims;
    output<<join(drugNameList,outputDelims)<<std::endl;
    int num = 0;
    for(const auto & protein : predictProteinsIndex){
      for(int j=0;j<drugNum;++j){
        tmpCalc.push_back(iterdrugSub2ProteinSub(j,protein));
      } // end of loop for j
      output<<existNameList[num]<<outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
      ++num;
    } // end of loop for protein
    return 0;
  } // end of functions.

  int EM::predictProteinsWithSubs(){
    std::cout<<"Now predict given proteins with subs against all the drugs, "
             <<"which are in our training data set." << std::endl;
    IntList predictProteinsIndex;
    std::vector<double> tmpCalc;
    nameList existNameList;
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as drug names.
    output<<"Names"<<outputDelims;
    output<<join(drugNameList,outputDelims)<<std::endl;
    double tmp;
    for(int i=0;i<predictProteinNameList_WithSubs.size();++i){
      for(int j=0;j<drugNum;++j){
        tmp = 0;
        for(auto const & m : (*drug2sub)[j]){
          for(auto const & n : predictProtein2SubList[i]){
            tmp += log(1 - (*drugSub2proteinSub)[m][n]);
          } // end of loop for n
        } // end of loop for m
        tmpCalc.push_back(exp(tmp));
      } // end of loop for j
      output<<predictProteinNameList_WithSubs[i]<<outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
    } // end of loop for protein
    return 0;
  } // end of functions.

  int EM::predictDrugsProteins(){
    std::cout<<"Now predict given drugs against given proteins." << std::endl;
    IntList predictDrugIndex;
    IntList predictProteinIndex;
    std::vector<double> tmpCalc;
    nameList existdrugNameList;
    nameList existproteinNameList;
    getIndexFromHash(drugName2Index, predictDrugNameList,
                     predictDrugIndex, existdrugNameList);
    getIndexFromHash(proteinName2Index, predictProteinNameList,
                     predictProteinIndex, existproteinNameList);
    if (existdrugNameList.size()<1) {
      std::cerr<<"No drug Index found, and quit." <<std::endl;
      return 1;
    } // end of if
    if (existproteinNameList.size()<1) {
      std::cerr<<"No protein Index found, and quit." << std::endl;
      return 1;
    } // end of if
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as protein names.
    output<<"Names"<<outputDelims;
    output<<join(existproteinNameList,outputDelims)<<std::endl;
    int num = 0;
    for(const auto & drug : predictDrugIndex){
      for(const auto & protein : predictProteinIndex){
        tmpCalc.push_back(iterdrugSub2ProteinSub(drug,protein));
      } // end of loop for protein
      output<<existdrugNameList[num] <<outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
      ++num;
    } // end of loop for drug
    return 0;
  } // end of functions.

  int EM::predictDrugsProteinsWithSubs(){
    std::cout<<"Now predict given drugs against given proteins with subs."
             << std::endl;
    IntList predictDrugIndex;
    std::vector<double> tmpCalc;
    nameList existdrugNameList;
    getIndexFromHash(drugName2Index, predictDrugNameList,
                     predictDrugIndex, existdrugNameList);
    if (existdrugNameList.size()<1) {
      std::cerr<<"No drug Index found, and quit." <<std::endl;
      return 1;
    } // end of if
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as protein names.
    output<<"Names"<<outputDelims;
    output<<join(predictProteinNameList_WithSubs,outputDelims)<<std::endl;
    int num = 0;
    double tmp;
    for(const auto & drug : predictDrugIndex){
      for(int j=0;j<predictProteinNameList_WithSubs.size();++j){
        tmp = 0;
        for(auto const & m : (*drug2sub)[drug]){
          for(auto const & n : predictProtein2SubList[j]){
            tmp += log(1 - (*drugSub2proteinSub)[m][n]);
          } // end of loop n
        } // end of loop for m
        tmpCalc.push_back(exp(tmp));
      } // end of loop for protein
      output<<existdrugNameList[num]<<outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
      ++num;
    } // end of loop for drug
    return 0;
  } // end of functions.

  int EM::predictDrugsWithSubsProteins(){
    std::cout<<"Now predict given drugs with subs against given proteins."
             << std::endl;
    IntList predictProteinIndex;
    std::vector<double> tmpCalc;
    nameList existproteinNameList;
    getIndexFromHash(proteinName2Index, predictProteinNameList,
                     predictProteinIndex, existproteinNameList);
    if (existproteinNameList.size()<1) {
      std::cerr<<"No protein Index found, and quit." <<std::endl;
      return 1;
    } // end of if
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as protein names.
    output<<"Names"<<outputDelims;
    output<<join(existproteinNameList,outputDelims)<<std::endl;
    double tmp;
    for(int i=0;i<predictDrugNameList_WithSubs.size();++i){
      for(auto const protein : predictProteinIndex){
        tmp = 0;
        for(auto const & m : predictDrug2SubList[i]){
          for(auto const & n : (*protein2sub)[protein]){
            tmp += log(1 - (*drugSub2proteinSub)[m][n]);
          } // end of loop n
        } // end of loop for m
        tmpCalc.push_back(exp(tmp));
      } // end of loop for protein
      output<<predictDrugNameList_WithSubs[i]<<outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
    } // end of loop for i
    return 0;
  } // end of functions.

  int EM::predictDrugsWithSubsProteinsWithSubs(){
    std::cout<<"Now predict given drugs with subs against given proteins with subs."
             << std::endl;
    std::vector<double> tmpCalc;
    // output the result.
    std::ofstream output (outPredictCPIsFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<outPredictCPIsFileName<<std::endl;
      return 1;
    } // end of if
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // print the first row as protein names.
    output<<"Names"<<outputDelims;
    output<<join(predictProteinNameList_WithSubs,outputDelims)<<std::endl;
    double tmp;
    for(int i=0;i<predictDrugNameList_WithSubs.size();++i){
      for(int j=0;j<predictProteinNameList_WithSubs.size();++j){
        tmp = 0;
        for(auto const & m : predictProtein2SubList[i]){
          for(auto const & n : predictProtein2SubList[j]){
            tmp += log(1 - (*drugSub2proteinSub)[m][n]);
          } // end of loop n
        } // end of loop for m
        tmpCalc.push_back(exp(tmp));
      } // end of loop for protein
      output<<predictDrugNameList_WithSubs[i]<<outputDelims;
      // transformed without static_cast should also work?
      output<<join(tmpCalc |
              transformed(static_cast<std::string(*)(double)>(std::to_string) ),
                   outputDelims)
            << std::endl;
      tmpCalc.clear();
    } // end of loop for i
    return 0;
  } // end of functions.

  int EM::outTrainResult(){
    writeMatrix(outDrugSub2ProteinSubFileName, *drugSub2proteinSub,
                outputDelims);
    return 0;
  } // end of function

  int EM::outTrainVariance(){
    writeMatrix(outVarDrugSub2proteinSubFileName, *vardrugSub2proteinSub,
                outputDelims);
    return 0;
  } // end of function

} // end of namespace gift.
