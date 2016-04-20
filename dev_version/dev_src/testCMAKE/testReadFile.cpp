#include<iostream>
#include<fstream>
#include<boost/program_options.hpp>

class testParameter{
public:
  testParameter(std::string);
  inline void outMember(){
    std::cout<<delims<<std::endl;
    std::cout<<double_fn<<std::endl;
    std::cout<<int_fn<<std::endl;
    std::cout<<task<<std::endl;
  }
private:
  std::string delims;
  double double_fn;
  int int_fn;
  std::string task;
};

testParameter::testParameter(std::string inputFile){
  std::cout<<"Now construct an object of testParameter..."<<std::endl;
  std::ifstream input;
  input.open(inputFile,std::ifstream::in);
  namespace po = boost::program_options;
  po::options_description desc("testParameter options");
  desc.add_options()
    ("delims",po::value<std::string>(&delims),"file sep")
    ("double_fn",po::value<double>(&double_fn),"double of fn")
    ("int_fn",po::value<int>(&int_fn),"int of fn")
    ("task",po::value<std::string>(&task),"task");
  po::variables_map vm;
  po::store(po::parse_config_file(input, desc),vm);
  po::notify(vm);
}

int main(int argc, char ** argv){
  namespace po = boost::program_options;
  po::options_description desc("Command Options");
  std::string inputFile;
  desc.add_options()
    ("config",po::value<std::string>(&inputFile),"inputFile name");
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  testParameter testCase(inputFile);
  testCase.outMember();
  return 0;
}
