// The source file for testing gift module.

// Author: Songpeng Zu
// Email: zusongpeng@gmail.com
// Date: 2016-03-06

// C/C++ system library
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
// Boost library
#define BOOST_REQUIRE_MODULE test_gift_functions
// When want to use self main.
//#define BOOST_REQUIRE_NO_MAIN
//#define BOOST_REQUIRE_ALTERNATIVE_INIT_API
#include<boost/test/included/unit_test.hpp>
// Third party library
// Own library
#include "gift.h" // Add gift.h to PATH.

namespace utf = boost::unit_test;

unsiged testStrEq(const std::string a, const std::string b){
  return a.compare(b);
}

// Init of class of parameters.
gift::parameters tmp_a; // test default inition.
delete tmp_a;

std::ifstream test_init ("test_init-predict.gift");
std::ifstream test_d2p_file ("test_drug2protein"); // QUESTION
std::ifstream test_p2d_file ("test_protein2domain"); // QUESTION
std::ifstream test_d2s_file ("test_drug2sub"); // QUESTION
std::ifstream test_s2d_file ("test_sub2domain"); // QUESTION

std::vector<std::vector<int> > test_d2p;
std::vector<std::vector<int> > test_p2d;
std::vector<std::vector<int> > test_d2s;
std::vector<std::vector<int> > test_s2d; // drugSub to proteinDomian.

gift::readMatrix<int>(test_d2p_file,test_d2p);
gift::readMatrix<int>(test_p2d_file,test_p2d);
gift::readMatrix<int>(test_d2s_file,test_d2s);
//gift::readMatrix<int>(test_s2d_file,test_s2d); // Task: prediction.

gift::parameters param (test_init);
test_init.close(); // end of file test_init-predict.gift

BOOST_AUTO_TEST_SUITE( test_gift_input_func)

//BOOST_AUTO_TEST_CASE( test_readMatrix ){
//  BOOST_REQUIRE( gift::readMatrix<int>(test_d2p_file,test_d2p) == 0);
//}
BOOST_AUTO_TEST_CASE( test_helpfunc ){
  BOOST_REQUIRE( gift::helpGift() == 0);
}
BOOST_AUTO_TEST_CASE( test_get_matrix_rowCol){
  gift::rowCol test_rowCol;
  BOOST_REQUIRE(gift::rowColFile("test_d2p_file",test_rowCol) == 0);
  BOOST_REQUIRE(test_rowCol.rowNum == 2);
  BOOST_REQUIRE(test_rowCol.colNum == 3);
}

BOOST_AUTO_TEST_SUITE_END() // end of test_gift_input_func

BOOST_AUTO_TEST_SUITE( test_class_parameters )

  BOOST_AUTO_TEST_CASE( test_dataMembers, *uft::tolerance(0.00001)){
  BOOST_REQUIRE(param.fn == 0.85);
  BOOST_REQUIRE(gift::testStrEq(gift::param.task,"train"));
  BOOST_REQUIRE(!param.loglikelyRecord);
  BOOST_REQUIRE(gift::testStrEq(gift::param.chemfpRec, "ComFP: PUBCHEM."));
  BOOST_REQUIRE(param.drugNum == 0);
}
BOOST_AUTO_TEST_CASE( test_funcMembers_setNumbers ){
  gift::rowCol param_set_rowCol;
  gift::rowColFile("test_drug2protein",param_set_rowCol)
  BOOST_REQUIRE(param.setDrugNum(param_set_rowCol.rowNum));
  BOOST_REQUIRE(param.setProteinNum(param_set_rowCol.colNum));

  gift::rowColFile("test_drug2sub",param_set_rowCol);
  BOOST_REQUIRE(param.setSubNum(param_set_rowCol.colNum));
  BOOST_REQUIRE(param.drugNum == param_set_rowCol.rowNum);

  gift::rowColFile("test_protein2domain",param_set_rowCol);
  BOOST_REQUIRE(param.setDomainNum(param_set_rowCol.colNum));
  BOOST_REQUIRE(param.proteinNum == param_set_rowCol.rowNum);
}

BOOST_AUTO_TEST_SUITE_END() // end of test_class_parameters.

BOOST_AUTO_TEST_SUITE( test_class_EM )

gift::EM test_em(param); // Global data ?

BOOST_AUTO_TEST_CASE( test_dataMembers ){
  BOOST_REQUIRE(test_em.iterationNum == param.iterationNum);
  BOOST_REQUIRE(test_em.drug2sub == NULL);
  test_em.setPointerDrug2Protein(test_d2p);
  BOOST_REQUIRE(test_em.drug2protein == &test_d2p);
  test_em.setPointerProtein2Sub(test_p2d);
  test_em.setPointerDrug2Sub(test_d2s);
  BOOST_REQUIRE(test_em.initEM());
}

BOOST_AUTO_TEST_CASE( test_memberFuncs ){
  BOOST_REQUIRE(EStep() == 0);
  BOOST_REQUIRE(MStep() == 0);
  BOOST_REQUIRE(trainEM() == 0);
  BOOST_REQUIRE(outTrain() == 0);
  BOOST_REQUIRE(loglikely() == 0);
// BOOST_REQUIRE( predictEM() == 0 );
}

BOOST_AUTO_TEST_SUITE_END() // end of test_class_EM

BOOST_AUTO_TEST_SUITE( test_output_func )

BOOST_AUTO_TEST_CASE( test_outRecord ){
  BOOST_REQUIRE(outRecord(param, test_em) == 0);
}

BOOST_AUTO_TEST_SUITE_END() // end of test_output_func
// When want to use self main.
//int main(int argc, char* argv[]){
//  some function...
//  return utf::unit_test_main(init_unit_test, argc, argv);
//}
