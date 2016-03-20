// functions without classes in the namespace gift.

// Include in order:
//   1. Matrix2Fingerprints;
//   2. readMatrix;
//   3. rowColFile;
//   4. writeMatrix;
//   5. helpGift;
//   6. outRecord;

// Libraries.
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

  int helpGift(){ // LACK OF DEFINITION
    return 0;
  } // end of function

  int outRecord(){ // LACK OF DEFINITION
    return 0;
  } // end of function

} // end of namespace gift.
