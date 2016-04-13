// The main file for the method of gift.

// Details can be found from the gift.h
// Author: Songpeng Zu
// Email: zusongpeng@gmail.com
// Date: 2016-04-05

// Load library.
#include<iostream>
#include<boost/program_options.hpp> // read the parameters
#include "gift.h"

namespace{
  const int SUCCESS = 0;
  const int ERROR_IN_COMMAND_LINE = 1;
  const int ERROR_UNHANDLED_EXCEPTION = 2;
  const int ERROR_IN_READFILE = 3;
  const int ERROR_IN_TASK = 4;
}

int main(int argc, char ** argv){
  if (argc < 2){
    gift::helpGift();
    return SUCCESS;
  } // end of if.
  // read program options.
  std::string configureFileName;
  // reference to: radmangames online post.
  try{
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
      ("help,h","Print help messages.")
      ("version,v","Print version information.")
      ("configure,c",po::value<std::string>(&configureFileName)->required(),
       "Read the configure file.");
    po::variables_map vm;
    try{
      po::store(po::parse_command_line(argc, argv, desc), vm);
      if(vm.count("help") || vm.count("-h")) {
        gift::helpGift();
      } // end of if
      if(vm.count("version") || vm.count("-v")) {
        std::cout<<"GIFT VERSION: "<<gift::version<<std::endl;
        std::cout<<"UPDATE TIME: "<<gift::updateTime<<std::endl;
      } // end of if

      po::notify(vm); // throw an error if there are any problems.

    } catch(po::error& e){
      std::cerr<<"ERROR: " <<e.what()<<std::endl<<std::endl;
      std::cerr<< desc <<std::endl;
      return ERROR_IN_COMMAND_LINE;
    } // end of catch
  } catch(std::exception& e){
    std::cerr<< "Unhandled  Exception reached the top of main: "
             <<e.what() <<", gift will now exit."<<std::endl;
    return ERROR_UNHANDLED_EXCEPTION;
  } // end of catch

  // run program based on task.
  try{
    gift::parameters getParameters(configureFileName);
    gift::EM EMgiftor(getParameters);
    if (getParameters.task.compare("predict")){
      // check is it enough?
      EMgiftor.predictEM();
      gift::outRecord(getParameters,EMgiftor);
    } else if (getParameters.task.compare("train")){
      // check is it enough?
      EMgiftor.trainEM();
      EMgiftor.varEM();
      // out train result.
      gift::outRecord(getParameters, EMgiftor);
      EMgiftor.outTrainResult();
      EMgiftor.outTrainVariance();
    } else {
      std::cout<<"NO TASK IS SPECIFIED."<<std::endl;
      return ERROR_IN_TASK;
    } // end of if else if else.
  } catch (const std::string & e){
    std::cerr<<"ERROR: "<<e<<std::endl;
    return ERROR_IN_READFILE;
  }
  return SUCCESS;
} // end of main
