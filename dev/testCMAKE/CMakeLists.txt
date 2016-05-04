CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT(TESTC++11)
SET(CMAKE_CXX_STANDARD 11)

# The two comments are not needed any more.
#SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
#SET(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++")

# SET(Boost_USE_STATIC_LIBS ON)
# SET(Boost_USE_MULTITHREADED ON)
# SET(Boost_USE_STATIC_RUNTIME OFF)
FIND_PACKAGE(Boost COMPONENTS thread program_options system exception)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "Boost library is " ${Boost_LIBRARIES})
#  LINK_DIRECTORIES(/usr/local/Cellar/boost/1.60.0_1/lib/)
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find boost library.")
ENDIF(Boost_FOUND)
#ADD_SUBDIRECTORY(src bin)
# test c++11 and boost library.
ADD_EXECUTABLE(testRead test_readMatrix.cpp)
# test programs_options reading from file.
#ADD_EXECUTABLE(testc11 testReadFile.cpp)
TARGET_LINK_LIBRARIES(testRead ${Boost_LIBRARIES})
