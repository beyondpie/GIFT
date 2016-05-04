// test read Matrix
// Songpeng Zu
// Tue May  3 22:02:12 CST 2016

#include<string>
#include<vector>
#include<iostream>
#include<fstream>
#include<boost/algorithm/string.hpp>
#include<boost/algorithm/string/join.hpp>
#include<boost/range/adaptor/transformed.hpp>

typedef std::vector<std::vector<double> > numericMatrix;

int printMatrix(const numericMatrix& fromMatrix){
  int rowNum = fromMatrix.size();
  int colNum = fromMatrix[0].size();
  for(int i=0;i<rowNum;++i){
    for(int j=0;j<colNum;++j){
      std::cout<<fromMatrix[i][j]<<",";
    } // end of loop for j
    std::cout<<std::endl;
  } // end of loop for i
  return 0;
} // end of function

int readMatrix(const std::string inputFile, numericMatrix& getMat,
               std::string delims){
  std::ifstream input (inputFile, std::ios::in);
  std::string line;
  std::vector<std::string> array;
  std::vector<double> tempRec;
  while(std::getline(input,line)){
    std::cout<<"Current Line is "<<line<<std::endl;
    boost::algorithm::split(array,line,boost::is_any_of(delims));
    int arraylen = array.size();
    //test
    std::cout<<"Current Arraylen in sub2sub is " <<arraylen<<std::endl;
    for(int i=0;i<arraylen;++i){
      std::string::size_type* idx = 0;
      tempRec.push_back(std::stod(array[i], idx));
    } // end of for
    getMat.push_back(tempRec);
    tempRec.clear();
  } // end of while
  return 0;
} // end of function

int main(){
  std::string filename("test_outdrugsub2proteinsubFile.pp");
  numericMatrix temp;
  std::string delims = ",";
  readMatrix(filename, temp, delims);
  printMatrix(temp);
  return 0;
} // end of main
