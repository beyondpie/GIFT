// Implementation of class parameter in namespace gift.

// Libraries
#include<boost/program_options.hpp>

#include "gift.h"

namespcae gift{
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
       "outFile prefix for predicted compound FPs-protein FPs interactions")
    po::variables_map vm;
    std::ifstream input;
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      // Allowed empty file with default values.
      input.open(inputFile, std::ifstream::in);
    //   if (input.peek() == std::ifstream::traits_type::eof()) {
    //     std::cerr<<inputFile<<" is empty." << std::endl;
    //     return 1;
    //   } // end of if.
      po::store(po::parse_config_file(input, desc), vm);
      po::notify(vm);
    } catch (std::ifstream::failure e) {
      std::cerr <<"Exceptions open/read file "<<configFile<<std::endl;
      return 1;
    } // end of try catch

    // default training data parameters.
    // they will be set when read data files.
    drugNum = 0;
    subNum = 0;
    domainNum = 0;
    proteinNum = 0;

    // notice the setting results.
    std::cout<<"parameters have been set."<<std::endl;
    return 0;
  } // end of class parameter constructor.
} // end of namespace gift.h
