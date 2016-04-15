CMAKE_MINIMUM_REQUIRED(VERSION 3.4)

SET(CMAKE_CXX_STANDARD 11)
MESSAGE(STATUS "The compiler is  ${CMAKE_CXX_COMPILER_ID}")
# Compiler-specific C++11 activation.
# IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
#   EXECUTE_PROCESS(
#     COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
#   IF(NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
#     MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
#   ENDIF()
# ELSEIF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
# ELSE ()
#   MESSAGE(FATAL_ERROR, "Your C++ compiler does not supprot c++11.")
# ENDIF()


PROJECT (GIFT CXX)
# The version number,
SET (GIFT_VERSION_MAJOR 2)
SET (GIFT_VERSION_MINOR 0)

# # configure a header file to pass some of the CMake settings to source code.
# configure_file (
#   "${PROJECT_SOURCE_DIR}/giftConfig.h.in"
#   "${PROJECT_BINARY_DIR}/giftConfig.h"
#   )

# # add the binary tree to the search path for including files.
# # so that can find giftConfig.h.
# INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

# SET(BOOST_ROOT /usr/local/include)
FIND_PACKAGE(Boost)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  INCLUDE_DIRECTORIES(/usr/local/include/boost)
  MESSAGE(STATUS "The boost library dir is " ${Boos_INCLUDE_DIRS})
  ADD_SUBDIRECTORY(src bin)


  SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

  MESSAGE(STATUS "gift binary dir is " ${PROJECT_BINARY_DIR})
  MESSAGE(STATUS "gift source dir is " ${PROJECT_SOURCE_DIR})

  TARGET_LINK_LIBRARIES(gift ${Boost_LIBRARIES})

ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Uable to find correct Boost library. Did you set BOOST_ROOT?")
ENDIF(Boost_FOUND)
