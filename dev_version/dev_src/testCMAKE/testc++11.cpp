#include<iostream>
#include<vector>
#include<boost/program_options.hpp>

int main(int argc, char ** argv){
  std::cout<<"hi,songpeng, test c++11 under clang with cmake."<<std::endl;
  std::vector<int> array(3,0);
  for(const auto & v : array){
    std::cout<<"output array element: "<<v<<std::endl;
  }
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h","hi, songpeng.");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if(vm.count("help") || vm.count("h")) {
    std::cout<<desc<<std::endl;
  }
  return 0;
}
