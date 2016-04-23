#include<iostream>
#include<vector>
#include<boost/program_options.hpp>
#include<boost/algorithm/string/join.hpp>

int x = 10;
std::vector<std::vector<double> > y;

int init_y(){
  std::vector<double> tmp;
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      tmp.push_back(rand());
    } // end of for
    y.push_back(tmp);
  } // end of for
  return 0;
} // end of function

class pointer_y{
public:
  int update_y();
  pointer_y() : pointer(&y){}
  int out_y();
private:
  std::vector<std::vector<double> > * pointer;
}; // end of class

int pointer_y::update_y(){
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      (*pointer).at(i).at(j) = i+j;
    } // end of loop
  } // end of loop
  return 0;
} // end of function

int pointer_y::out_y(){
  std::cout<<"Print intMatrix y..."<<std::endl;
  for(int i =0;i<3;++i){
    for(int j=0;j<3;++j){
      std::cout<<(*pointer)[i][j]<<std::endl;
    } // end of loop
  } // end of loop
  return 0;
} // end of function


void changeX(){
  x = 5;
} // end of function

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
  std::vector<std::string> list;
  list.push_back("algorithm");
  list.push_back("stringjoin");
  std::string joined = boost::algorithm::join(list,";");
  std::cout<<joined<<std::endl;

  std::cout<<"Init x is "<<x<<std::endl;
  changeX();
  std::cout<<"Now x is "<<x<<std::endl;

  init_y();
  pointer_y a_point_y;
  a_point_y.out_y();
  a_point_y.update_y();
  a_point_y.out_y();
  return 0;
}
