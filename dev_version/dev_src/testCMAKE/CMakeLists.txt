CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT(TESTC++11)
SET(CMAKE_CXX_STANDARD 11)

# The two comments are not needed any more.
#SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
#SET(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++")

SET(Boost_USE_STATIC_LIBS ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
FIND_PACKAGE(Boost COMPONENTS program_options)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "Boost library is " ${Boost_LIBRARIES})
#  LINK_DIRECTORIES(/usr/local/Cellar/boost/1.60.0_1/lib/)
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find boost library.")
ENDIF(Boost_FOUND)

#INCLUDE_DIRECTORIES(/usr/local/Cellar/boost/1.60.0_1/include)
ADD_EXECUTABLE(testc11 testc++11.cpp)
#TARGET_LINK_LIBRARIES(testc11 /usr/local/Cellar/boost/1.60.0_1/lib)
TARGET_LINK_LIBRARIES(testc11 ${Boost_LIBRARIES})
