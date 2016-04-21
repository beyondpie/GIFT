// Implementation of class parameter in namespace gift.

// Libraries
#include "parameterGift.hpp"

namespace gift{
  inline const char* BoolToString (bool b) {
    return b ? "true" : "false";
  } // end of function BoolToString
  parameters::parameters(const std::string configFile) throw(std::string) {
    std::cout<<"Now set parameters with configFile."<<std::endl;
    // use boost program_options to read configs from a given file.
    namespace po = boost::program_options;
    po::options_description desc("GIFT Parameter options");
    desc.add_options()
      // input data file name
      ("drug2proteinFileName",
       po::value<std::string>(&drug2proteinFileName),
       "file name for drug protein interactions")
      ("drug2subFileName",po::value<std::string>(&drug2subFileName),
       "file name for drug to substructure")
      ("protein2subFileName",po::value<std::string>(&protein2subFileName),
       "file name for protein to substructure")
      ("drugSub2proteinSubFileName",
       po::value<std::string>(&drugSub2proteinSubFileName),
       "file name for drugSub to proteinSub interaction probability")
      // input name list file name
      ("drugNameListFile",po::value<std::string>(&drugNameListFile),
       "file name for drug names")
      ("drugSubNameListFile",po::value<std::string>(&drugSubNameListFile),
       "file name for drug substructures names")
      ("proteinNameListFile",po::value<std::string>(&proteinNameListFile),
       ("file name for protein names"))
      ("proteinSubNameListFile",po::value<std::string>(&proteinSubNameListFile),
       "file name for protein substructures names")
      // input parameters for EM
      ("alphaEB",po::value<double>(&alphaEB)->default_value(0.05),
       "parameter for Empirical Bayesian estimates for initEM")
      ("betaEB",po::value<double>(&betaEB)->default_value(0.05),
       "parameter for Empirical Bayesian estimates for initEM")
      ("fp", po::value<double>(&fp)->default_value(0.85),
       "false positive rate")
      ("fn", po::value<double>(&fn)->default_value(0.0001),
       "false negative rate")
      ("threadNum", po::value<int>(&thread)->default_value(1),
       "thread number for EM")
      ("EMIterationNum", po::value<int>(&iterNum)->default_value(300),
       "iteration numbers for EM")
      ("task", po::value<std::string>(&task)->default_value("train"),
       "run gift for train or predict")
      ("loglikelyRecord",
       po::value<bool>(&loglikelyRecord)->default_value(false),
       "whether or not to record the loglikely in every step")
      ("inputDelims",po::value<std::string>(&inputDelims)->default_value("\t,"),
       "sep character for input files")
      // input file version information.
      ("chemFingerPrintRecord",
       po::value<std::string>(&chemfpRec)->default_value("ComFP: PUBCHEM"),
       "source and version of chemical fingerprints")
      ("proteinFingerPrintRecord",
       po::value<std::string>(&proteinfpRec)->default_value("Pfam: 2011-07"),
       "source and version of protein fingerprints/domains")
      ("comProteinInteractionRecord",
       po::value<std::string>(&CPIsRec)->default_value("DrugBank: 2011-07"),
       "source and version of compound-protien interactions")
      // input file names for prediction
      ("predictDrugsFileName",po::value<std::string>(&predictDrugsFileName),
       "file name for drug names used for prediciton by gift")
      ("predictProteinFileName",po::value<std::string>(&predictProteinsFileName),
       "file name for protein names used for prediction by gift")
      ("predictDrugsFileName_WithSubs",
       po::value<std::string>(&predictDrugsFileName_WithSubs),
       "file name for drugs names together with their substructures.")
      ("predictProteinsFileName_WithSubs",
       po::value<std::string>(&predictProteinsFileName_WithSubs),
       "file name for protein names together with their substructures.")
      // output file name and format
      ("outputDelims",po::value<std::string>(&outputDelims)->default_value("\t,"),
       "sep character for output files")
      ("outRecordFileName",
       po::value<std::string>(&outRecordFileName)->default_value("CPIs"),
       "file name for output records")
      ("outPredictCPIsFileName",
       po::value<std::string>(&outPredictCPIsFileName),
       "file name for output CPIs")
      ("outDrugSub2ProteinSubFileName",
       po::value<std::string>(&outDrugSub2ProteinSubFileName),
       "file name for output drugSub2proteinSub")
      ("outVarDrugSub2proteinSubFileName",
       po::value<std::string>(&outVarDrugSub2proteinSubFileName),
       "file name for output variance of drugSub2proteinSub");
    po::variables_map vm;
    std::ifstream input;
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(configFile, std::ifstream::in);
      if (input.peek() == std::ifstream::traits_type::eof()) {
        std::cerr<<configFile<<" is empty." << std::endl;
      } // end of if.
      po::store(po::parse_config_file(input, desc), vm);
      po::notify(vm);
    } catch (std::ifstream::failure e) {
      std::cerr <<"Exceptions open/read file "<<configFile<<std::endl;
    } // end of try catch

    // load global data for gift.
    Matrix2Fingerpints(drug2proteinFileName,drug2proteinList,inputDelims);
    Matrix2Fingerpints(protein2subFileName,protein2domainList,inputDelims);
    Matrix2Fingerpints(drug2subFileName,drug2subList,inputDelims);
    InitDrugSub2ProteinSub();
    InitVarDrugSub2ProteinSub();
    // init sub2drugList and domain2proteinList.
    IntList tmpIntArray;
    for(int i=0;i<subNum;++i){sub2drugList.push_back(tmpIntArray);}
    for(int i=0;i<domainNum;++i){domain2proteinList.push_back(tmpIntArray);}
    Matrix2FingerprintsByColumn(drug2subFileName,sub2drugList,subNum,inputDelims);
    Matrix2FingerprintsByColumn(protein2subFileName,domain2proteinList,
                                domainNum, inputDelims);

    // load NameList.
    InitDrugName2Index();
    InitProteinName2Index();
    InitDrugSubNameList();
    InitProteinSubNameList();

    // load predicted name list and possible subs if task is prediction.
    InitPredictParameters(); // throw std::string.

    // default training data parameters.
    // they will be set when read data files.
    rowCol tmp;
    rowColFile(drug2subFileName,tmp,inputDelims);
    drugNum = tmp.rowNum;
    subNum = tmp.colNum;
    rowColFile(protein2subFileName,tmp,inputDelims);
    domainNum = tmp.colNum;
    proteinNum = tmp.rowNum;

    // print the setting results.
    std::cout<<"parameters have been set."<<std::endl;
    for (const auto& it : vm){
      std::cout<< it.first.c_str() << ": ";
      auto& value = it.second.value(); // return boost::any reference type.
      // any_cast use the any * as input and return the pointer with type infor.
      if (auto v = boost::any_cast<int>(&value) ) {
        std::cout<< *v <<std::endl;
      } else if (auto v = boost::any_cast<double>(&value) ) {
        std::cout<< *v <<std::endl;
      } else if (auto v = boost::any_cast<bool>(&value) ) {
        std::cout<< BoolToString(*v) <<std::endl;
      } else if (auto v = boost::any_cast<std::string>(&value) ) {
        std::cout<< *v <<std::endl;
      } else {
        std::cout<< "Error type"<<std::endl;
      } // end of if
    } // end of for
    std::cout<<"drugNum: " <<drugNum<<std::endl;
    std::cout<<"subNum: "<<subNum<<std::endl;
    std::cout<<"domainNum: "<<domainNum<<std::endl;
    std::cout<<"proteinNum: "<<proteinNum<<std::endl;
  } // end of class parameter constructor.

  int parameters::InitDrugSub2ProteinSub(){
    // This function must be run after class parameter initionlization.
    std::cout<< "Initialize the drugSub2proteinSub Matrix." << std::endl;
    if (task.compare("predict")) {
      readMatrix(drugSub2proteinSubFileName,drugSub2proteinSubMatrix,
                 inputDelims);
      std::cout<< "Finish: read from file."<<std::endl;
    } else {
      std::vector<double> assoTmp; // temp array based on the assocaiton method.
      int N = 0;
      int subNumTmp = 0;
      int I = 0;
      std::vector<int>::iterator it;
      for (int i=0;i<subNum;++i){
        subNumTmp = sub2drugList[i].size();
        for (int j=0;j<domainNum;++j){
          N = domain2proteinList[j].size() * subNumTmp;
          for (const auto drug : sub2drugList[i]){
            for (const auto protein : domain2proteinList[j]){
              it = std::find(drug2proteinList[drug].begin(),
                             drug2proteinList[drug].end(), protein);
              I += it==drug2proteinList[drug].end() ? 0 : 1;
            } // end of loop protein
          } // end of loop drug
          // revise association method with Emiprical Bayes.
          assoTmp.push_back((I+alphaEB)/(alphaEB+betaEB+N));
        } // end of loop j
        drugSub2proteinSubMatrix.push_back(assoTmp);
        assoTmp.empty();
      } // end of loop i
      std::cout<<"Finish: initialize with associatiom method and emprical Bayes."
               <<std::endl;
    } // end of if else
    return 0;
  } // end of function

  int parameters::InitVarDrugSub2ProteinSub(){
    if (task.compare("predict")) {
      std::cout<< "Job is to predict, skip init variance matrix." << std::endl;
    } else {
      std::cout<<"Initialize the variance matrix for drugSub2ProteinSub."
               <<std::endl;
      std::vector<double> tmpArray (domainNum,0.0);
      for(int i=0;i<subNum;++i){
        vardrugSub2proteinSubMatrix.push_back(tmpArray);
      } // end of loop for i.
    } // end of if else
    return 0;
  } // end of function

  int parameters::InitDrugName2Index() {
    readNameListFromFile(drugNameListFile,drugNameList);
    readName2IndexHash(drugNameList,drugName2Index);
    return 0;
  } // end of function

  int parameters::InitProteinName2Index() {
    readNameListFromFile(proteinNameListFile, proteinNameList);
    readName2IndexHash(proteinNameList,proteinName2Index);
    return 0;
  } // end of function

  int parameters::InitDrugSubNameList() {
    readNameListFromFile(drugSubNameListFile,drugSubNameList);
    return 0;
  } // end of function

  int parameters::InitProteinSubNameList() {
    readNameListFromFile(proteinSubNameListFile,proteinSubNameList);
    return 0;
  } // end of function

  int parameters::InitPredictParameters() throw(std::string) {
    // when task is predict, we use this function to init the corresponding
    // parameters.
    bool checkStatus = false;
    if (!predictDrugsFileName.empty()){
      checkStatus = true;
      readNameListFromFile(predictDrugsFileName,predictDrugNameList);
    } // end of if for drugfile
    if (!predictProteinNameList.empty()){
      checkStatus = true;
      readNameListFromFile(predictDrugsFileName,predictProteinNameList);
    } // end of if for proteinfile
    if (!predictDrugNameList_WithSubs.empty()){
      checkStatus  = true;
      // NEED FUNCTION.
      readNameMatrixFromFile(predictDrugsFileName_WithSubs,
                             predictDrugNameList_WithSubs,
                             predictDrug2SubList);
    } // end of if for drug_withsubs file.
    if (!predictProteinNameList_WithSubs.empty()){
      checkStatus = true;
      // NEED FUNCTION.
      readNameMatrixFromFile(predictProteinsFileName_WithSubs,
                             predictProteinNameList_WithSubs,
                             predictProtein2SubList);
    } // end of if for protein_withsubs file.
    if (!checkStatus){
      throw("Task for prediction, but no predicted files!");
    } // end of if for checkStatus.
    return 0;
  } // end of function

} // end of namespace gift.h
