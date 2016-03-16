// Libraries.
#include<boost/algorihtm/string.hpp>
#include "gift.h"

namespace gift{

  template <typename T> int readMatrix(const std::ifstream& inputMatrix,
                                       std::vector<std::vector<T> >& readM,
                                       std::string delims = "\t,"){
    std::string line;
    std::vector<T> tmp_array;
  }

  int rowColFile(const std::string inputFile, rowCol& matrixRec,
                 std::string delims="\t,"){
    std::ifstream input;
    std::string line;
    int count = 0;
    std::vector<string> array;

    input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      input.open(inputFile,std::ifstream::in);

      if (input.peek() == std::ifstream::traits_type::eof()){
        std::cerr <<inputFile<<" is empty."<<endl;
        return 1;
      }

      std::getline(input,line);
      ++count;
      boost::split(array, line, boost::is_any_of(delims));
      matrixRec.colNum = array.size();
      // string getline func over istream.
      while(std::getline(input,line)){
        // QUESTION: how about empty line?
        ++count;
      } // end of while
      input.close();

      matrixRec.rowNum = count;
      matrixRec.colNum = array.size();
    } //end of try
    catch (std::ifstream::failure e) {
      std::cerr << "Exceptions open/read/close file "<<inputFile<<std::endl;
      return 1;
    } // end of catch
    return 0;
  } // end of function rowColFile.

}
