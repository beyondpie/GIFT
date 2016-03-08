// The source file for testing gift module.

// Author: Songpeng Zu
// Email: zusongpeng@gmail.com
// Date: 2016-03-06

// C/C++ system library
#include<iostream>
#include<fstream>
#include<string>
// Boost library
#define BOOST_TEST_MODULE test_gift_functions
// When want to use self main.
//#define BOOST_TEST_NO_MAIN
//#define BOOST_TEST_ALTERNATIVE_INIT_API
#include<boost/test/included/unit_test.hpp>
// Third party library
// Own library
#include "gift.h" // Add gift.h to PATH.

namespace utf = boost::unit_test;

unsiged testStrEq(const std::string a, const std::string b){
  return a.compare(b);
}

// Init of class of parameters.
gift::parameters tmp_a;
std::ifstream test_init ("test_init-predict.gift");
gift::parameters param (test_init);
test_init.close();

BOOST_AUTO_TEST_SUITE( test_class_parameters )
BOOST_AUTO_TEST_CASE( test_dataMembers, *uft::tolerance(0.00001)){
  BOOST_TEST(param.fn == 0.85);
  BOOST_TEST(testStrEq(param.task,"train"));
  BOOST_TEST(!param.loglikelyRecord);
  BOOST_TEST(testStrEq(param.chemfpRec, "ComFP: PUBCHEM."));
  BOOST_TEST(param.drugNum == 0);
}
BOOST_AUTO_TEST_CASE( test_funcMembers ){
}
BOOST_AUTO_TEST_SUITE_END() // end of test_class_parameters.

BOOST_AUTO_TEST_SUITE( test_class_EM )
BOOST_AUTO_TEST_CASE( test_dataMembers ){
}
BOOST_AUTO_TEST_SUITE_END()

// When want to use self main.
//int main(int argc, char* argv[]){
//  some function...
//  return utf::unit_test_main(init_unit_test, argc, argv);
//}
