CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

SET(CMAKE_CXX_STANDARD 11)

# The two comments are not needed any more.
#SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
#SET(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++")

PROJECT (TESTC++11)
ADD_EXECUTABLE(testc11 testc++11.cpp)
