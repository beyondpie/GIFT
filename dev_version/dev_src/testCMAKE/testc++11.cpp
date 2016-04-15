#include<iostream>
#include<vector>

int main(){
  std::cout<<"hi,songpeng, test c++11 under clang with cmake."<<std::endl;
  std::vector<int> array(10,0);
  for(const auto & v : array){
    std::cout<<"output array element: "<<v<<std::endl;
  }
  return 0;
}
