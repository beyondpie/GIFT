// Implementation of class parameter in namespace gift.

// Libraries
#include<boost/program_options.hpp>
#include<boost/any.hpp>

#include "gift.h"

namespace gift{
  inline const char* BoolToString (bool b) {
    return b ? "true" : "false";
  } // end of function BoolToString

  parameters::parameters(const std::string configFile){
    std::cout<<"Now set parameters with configFile."<<std::endl;
    // use boost program_options to read configs from a given file.
    namespace po = boost::program_options;
    po::options_description desc("GIFT Parameter options");
    // Does it work for class members?
    desc.add_options()
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
      ("chemFingerPrintRecord",
       po::value<std::string>(&chemfpRec)->default_value("ComFP: PUBCHEM"),
       "source and version of chemical fingerprints")
      ("proteinFingerPrintRecprd",
       po::value<std::string>(&proteinfpRec)->default_value("Pfam: 2011-07"),
       "source and version of protein fingerprints/domains")
      ("comProteinInteractionRecord",
       po::value<std::string>(&CPIsRec)->default_value("DrugBank: 2011-07"),
       "source and version of compound-protien interactions")
      ("outFilePrefixCPIs",
       po::value<std::string>(&outFilePrefixCPIs)->default_value("CPIs"),
       "outFile prefix for predicted drug-protein interactions")
      ("outFilePrefixTrain",
       po::value<std::string>(&outFilePrefixTrain)->default_value("chemFP2proFP"),
       "outFile prefix for predicted compound FPs-protein FPs interactions");
    po::variables_map vm;
    std::ifstream input;
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      // Allowed file's content is empty.
      input.open(configFile, std::ifstream::in);
    //   if (input.peek() == std::ifstream::traits_type::eof()) {
    //     std::cerr<<inputFile<<" is empty." << std::endl;
    //     return 1;
    //   } // end of if.
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
      //
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
} // end of namespace gift.h
