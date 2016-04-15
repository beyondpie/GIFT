// functions without classes in the namespace gift.

// Libraries.
#include<ctime> // for record timing.
#include<chrono> // for record time.
#include<boost/algorithm/string.hpp>
#include<boost/algorithm/string/join.hpp>
#include<boost/range/adaptor/transformed.hpp>

#include "gift.h"

namespace gift{

  int Matrix2Fingerpints(const std::string inputFile, IntArrayList & getFp,
                         std::string delims){
    std::ifstream input;
    std::string line;
    std::vector<std::string> array;
    std::vector<int> tempRec;

    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(inputFile, std::ifstream::in);
      if (input.peek() == std::ifstream::traits_type::eof()){
        std::cerr <<inputFile <<" is empty. " <<std::endl;
        return 1;
      } // end of if
      while (std::getline(input,line)) {
        boost::algorithm::split(array,line,boost::is_any_of(delims));
        int arraylen = array.size();
        for (int i=0;i<arraylen;++i) {
          if (array[i].compare("1") == 0) {
            tempRec.push_back(i);
          }// end of if
        } // end of for
        getFp.push_back(tempRec);
        tempRec.clear();
      } // end of while
    } catch (std::ifstream::failure e) {
      std::cerr <<"Exceptions open/read file "<<inputFile<<std::endl;
      return 1;
    } // end of catch
    return 0;
  } // end of function.

  int readMatrix(const std::string inputFile, numericMatrix& getMat,
                 std::string delims){
    std::ifstream input;
    std::string line;
    std::vector<std::string> array;
    std::vector<double> tempRec;

    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(inputFile, std::ifstream::in);
      if (input.peek() == std::ifstream::traits_type::eof()){
        std::cerr <<inputFile <<" is empty. " <<std::endl;
      } // end of if
      while(std::getline(input,line)){
        boost::algorithm::split(array,line,boost::is_any_of(delims));
        int arraylen = array.size();
        for(int i=0;i<arraylen;++i){
          std::string::size_type* idx = 0;
          tempRec.push_back(std::stod(array[i], idx));
        } // end of for
        getMat.push_back(tempRec);
        tempRec.clear();
      } // end of while
    } catch (std::ifstream::failure e) {
      std::cerr <<"Exceptions open/read file "<<inputFile<<std::endl;
      return 1;
    } // end of catch
    return 0;
  }// end of function.

  int rowColFile(const std::string inputFile, rowCol& matrixRec,
                 std::string delims){
    std::ifstream input;
    std::string line;
    int count = 0;
    std::vector<std::string> array;

    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(inputFile,std::ifstream::in);
      if (input.peek() == std::ifstream::traits_type::eof()){
        std::cerr <<inputFile<<" is empty."<<std::endl;
        return 1;
      }

      std::getline(input,line);
      ++count;
      boost::algorithm::split(array, line, boost::is_any_of(delims));
      matrixRec.colNum = array.size();
      // string getline func over istream.
      while(std::getline(input,line)){
        // QUESTION: how about empty line?
        ++count;
      } // end of while
      input.close();

      matrixRec.rowNum = count;
      matrixRec.colNum = array.size();
    } catch (std::ifstream::failure e) {
      std::cerr << "Exceptions open/read file "<<inputFile<<std::endl;
      return 1;
    } // end of catch
    return 0;
  } // end of function.

  int writeMatrix(const std::string outFileName, numericMatrix& resultMat,
                  std::string delims){
    std::ofstream output (outFileName,std::ofstream::out);
    if (output.is_open()) {
      using boost::algorithm::join;
      using boost::adaptors::transformed;
      for (numericMatrix::iterator it = resultMat.begin();
           it != resultMat.end(); ++it) {
        output << join(*it | transformed(static_cast<std::string(*)(double)>
                                         (std::to_string) ), delims);
      }// end of for
      output.close();
    } else {
      std::cerr<< "Error opening file " <<outFileName<<std::endl;
      return 1;
    } // end of if else
    return 0;
  }// end of function

  int readNameListFromFile(const std::string inputFile, nameList& tonameList){
    // each line in the file represents one name.
    // line should end with "\n", not "[\r\t]\n"
    std::ifstream input;
    std::string line;
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(inputFile,std::ifstream::in);
      if (input.peek() == std::ifstream::traits_type::eof() ){
        std::cerr << inputFile <<" is empty. "<<std::endl;
        return 1;
      } // end of if
      while (std::getline(input,line)) {
        tonameList.push_back(line);
      } // end of while
      input.close();
    } catch (std::ifstream::failure e) {
      std::cerr<<"Exceptions open/read file " << inputFile<<std::endl;
    } // end of try catch
    return 0;
  } // end of function

  int readNameMatrixFromFile(const std::string inputFile, nameList& tonameList,
                             IntArrayList& getFP, std::string delims){
    std::ifstream input;
    std::string line;
    std::vector<std::string> array;
    std::vector<int> tempRec;

    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(inputFile, std::ifstream::in);
      if (input.peek() == std::ifstream::traits_type::eof()) {
        std::cerr<<inputFile<<" is empty. "<<std::endl;
        return 1;
      } // end of if
      while(std::getline(input,line)) {
        boost::algorithm::split(array,line,boost::is_any_of(delims));
        tonameList.push_back(array[0]); // first column is name.
        int arraylen = array.size();
        for(int i=1;i<arraylen;++i){
          if(array[i].compare("1") == 0) {
            tempRec.push_back(i-1); // Use i-1, since first column is name.
          } // end of if
        } // end of loop for i.
        getFP.push_back(tempRec);
        tempRec.clear();
      } // end of while
    } catch (std::ifstream::failure e) {
      std::cerr<<"Exceptions open/read file "<<inputFile<<std::endl;
      return 1;
    } // end of catch
    return 0;
  } // end of function

  int readName2IndexHash(const nameList fromNameList,
                         name2IndexHash& name2Index){
    // fromNameList should be in order.
    if (fromNameList.empty()){
      std::cerr<<"The fromNameList is empty. "<<std::endl;
      return 1;
    } // end of if
    int recordIndex = 0;
    for(const auto fromName : fromNameList){
      name2Index.insert(std::pair<std::string,int>(fromName,recordIndex));
      ++recordIndex;
    } // end of loop fromNameList
    return 0;
  } // end of function

  int getIndexFromHash(const name2IndexHash& name2Index,
                       const nameList fromNameList,
                       IntList & toIndexList,
                       nameList & existNameList){
    for(const auto fromName : fromNameList){
      if (name2Index.find(fromName) != name2Index.end()){
        existNameList.push_back(fromName);
        toIndexList.push_back( (name2Index.find(fromName))->second);
      } else {
        std::cout<<"Cannot find the key " << fromName
                 <<" fromNameList. Continue..." <<std::endl;
      } // end of if else
    } // end of loop fromNameList
    return 0;
  } // end of function

  int helpGift(){
    // Output gift information and useness to standard output.
    // Basic information about gift.
    std::cout<<"Gift is used to predict compound-protein interactions based on "
             <<std::endl;
    std::cout<<"their substructures interactions."<<std::endl;
    std::cout<<"It is also used to infer the substructres interactions from"
             <<std::endl;
    std::cout<<" the known drug-protein interactions."<<std::endl;
    std::cout<<"If you want to know more about gift, please read the paper: "
             <<std::endl;
    std::cout<<"Global Optimization-based Inference of Chemogenomic Features "
             <<std::endl;
    std::cout<<"from Drug-Target Interactions, which is published  "<<std::endl;
    std::cout<<"on Bioinformatics, 2015. "<<std::endl;
    std::cout<<"Author: "<<author<<std::endl;
    std::cout<<"Email: " <<email<<std::endl;
    std::cout<<"Current version: "<<version<<std::endl;
    std::cout<<"Last update time: "<<updateTime<<std::endl;
    std::cout<<"You can get the C++ source code from: "<<std::endl;
    std::cout<<"https://github.com/songpeng/GIFT" << std::endl;
    std::cout<<std::endl;

    //Input parameters
    std::cout<<"--help | -h to show the help information of gift."<<std::endl;
    std::cout<<"--version | -v to show the version information of gift."
             <<std::endl;
    std::cout<<"Gift need one configure file for its running."<<std::endl;
    std::cout<<"Please use --config to tell gift the configure file name."
             <<std::endl;
    std::cout<<"The content in the configure file are listed below: "<<std::endl;

    // configure file information.
    std::cout<<"[INPUT DATA FILE NAMES]" <<std::endl;
    std::cout<<"drug2proteinFileName=<string> : "
             <<"file name for drug protein interactions" <<std::endl;
    std::cout<<"drug2subFilename=<string> : "
             <<"file name for drug to substructure" <<std::endl;
    std::cout<<"protein2subFileName=<string> : "
             <<"file name for protein to substructure" << std::endl;
    std::cout<<"drugSub2proteinSubfilename=<string> : "
             <<"file name for drugSub to proteinSub interaction probability."
             <<std::endl;
    std::cout<<"drugNameListFile=<string> : "
             <<"file name for drug names" << std::endl;
    std::cout<<"drugSubNameListFile=<string> : "
             <<"file name for drug substructures names." << std::endl;
    std::cout<<"proteinNameListFile=<string> : "
             <<"file name for protein names" << std::endl;
    std::cout<<"proteinSubNameListFile=<string> : "
             <<"file name for protein substructures names." << std::endl;

    std::cout<<"[INPUT PARAMETERS FOR EM ALGORITHM]" <<std::endl;
    std::cout<<"alphaEB=<double> : "
             <<"parameter for Empricial Bayesian estimates for initEM."
             <<std::endl;
    std::cout<<"betaEB=<double> : "
             <<"parameter for Empricial Bayesian estimates for initEM."
             <<std::endl;
    std::cout<<"fp=<double> : "<<"false positive rate"<<std::endl;
    std::cout<<"fn=<double> : "<<"false negative rate"<<std::endl;
    std::cout<<"threadNum=<int> : "<<"thread number for EM." <<std::endl;
    std::cout<<"EMIterationNum=<int> : "<<"iteration numbers/steps for EM."<<std::endl;
    std::cout<<"task=<string> : "<<"run gift for [train] or [predict]."<<std::endl;
    std::cout<<"loglikelyRecord=<string> : "<<
      "record [true] or not [false] the loglikely in in every step."<<std::endl;
    std::cout<<"inputDelims=<string> : "
             <<"sep character for input files, such as  '\t',',' "<<std::endl;

    std::cout<<"[INPUT FILE VERSION INFORMATION]"<<std::endl;
    std::cout<<"chemFingerPrintRecord=<string> : "
             <<"source and version of chemical fingerprints."<<std::endl;
    std::cout<<"proteinFingerPrintRecord=<string> : "
             <<"source and version of protein fingerpints/domains."<<std::endl;
    std::cout<<"comProteinInteractionRecord=<string> : "
             <<"source and version of compound-protein interactions."<<std::endl;

    std::cout<<"[INPUT FILE NAME FOR PREDICTION]"<<std::endl;
    std::cout<<"predictDrugsFileName=<string> : "
             <<"file name for drug names used for prediction by gift."<<std::endl;
    std::cout<<"predictProteinsFileName=<string> : "
             <<"file name for protein names used for prediction by gift."<<std::endl;
    std::cout<<"predictDrugsFileName_WithSubs=<string> : "
             <<"file name for drug names together with their substructures."
             <<std::endl;
    std::cout<<"predictProteinsFileName_WithSubs=<string> : "
             <<"file name for protein names together with their substructures."
             <<std::endl;

    std::cout<<"[OUTPUT FILE NAME AND FORMAT]"<<std::endl;
    std::cout<<"outputDelims=<string> : "
             <<"sep character for output files." <<std::endl;
    std::cout<<"outRecordFileName=<string> : "<<"file name for output records."
             <<std::endl;
    std::cout<<"outPredictCPIsFileName=<string> : "
             <<"file name for output CPIs." <<std::endl;
    std::cout<<"outDrugSub2ProteinSubFileName=<string> : "
             <<"file name for output drugSub2proteinSub matrix."<<std::endl;
    std::cout<<"outVarDrugSub2proteinSubFileName=<string> : "
             <<"file name for output variance of drugSub2proteinSub." <<std::endl;

    return 0;
  } // end of function

  int outRecord(parameters & EMparameters, EM& EMgift){
    std::ofstream output (EMparameters.outRecordFileName,std::ofstream::out);
    if (!output.is_open()){
      std::cerr<<"Error open file "<<EMparameters.outRecordFileName<<std::endl;
      return 1;
    }// end of if
    // Basic information about gift.
    output<<"The author of gift is  " << author <<std::endl;
    output<<"Contact information: " << email <<std::endl;
    output<<"Current gift's version is "<< version <<std::endl;
    output<<"Update time is " << updateTime <<std::endl;
    // Running information.
    std::chrono::system_clock::time_point timePos =
      std::chrono::system_clock::now();
    std::time_t timePosT = std::chrono::system_clock::to_time_t(timePos);
    output<<"The job destination is "<<EMparameters.task;
    output<<", which is finished at "<<std::ctime(&timePosT) <<std::endl;
    output.close();
    return 0;
  } // end of function

  int Matrix2FingerprintsByColumn(const std::string inputFile,
                                  IntArrayList& getFP, int rowNum,
                                  std::string delims){
    std::ifstream input;
    std::string line;
    std::vector<std::string> array;
    int linenum = 0;
    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(inputFile,std::ifstream::in);
      if(input.peek() == std::ifstream::traits_type::eof()){
        std::cerr<< inputFile << " is empty. "<<std::endl;
        return 1;
      } // end of if

      // init getFP first.
      for(int i=0;i<rowNum;++i) {
        std::vector<int> tmpArray;
        getFP.push_back(tmpArray);
      } // end of loop

      // read file.
      while (std::getline(input,line)){
        boost::algorithm::split(array,line,boost::is_any_of(delims));
        int arraylen = array.size();
        for (int i=0;i<arraylen;++i){
          if (array[i].compare("1") == 0){
            getFP[i].push_back(linenum);
          } // end of if
        } // end of loop for i
        linenum += 1;
      } // end of while for file read.
    } catch (std::ifstream::failure e) {
      std::cerr<<"Exceptions open/read file "<<inputFile<<std::endl;
      return 1;
    } // end of catch
    return 0;
  } // end of function

} // end of namespace gift.
