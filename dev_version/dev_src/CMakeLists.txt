CMAKE_MINIMUM_REQUIRED(VERSION 3.4)

PROJECT (GIFT)
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
FIND_PACLAGE(Boost)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "The boost library dir is " ${Boost_LIBRARIES})
  LINK_DIRECTORIES(/usr/local/include) # ??
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Uable to find correct Boost library. Did you set BOOST_ROOT?")
ENDIF(Boost_FOUND)
# set source file list.
# SET(SRC_LIST main.cpp)
ADD_SUBDIRECTORY(src bin)

# set target dir.
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

MESSAGE(STATUS "gift binary dir is " ${PROJECT_BINARY_DIR})
MESSAGE(STATUS "gift source dir is " ${PROJECT_SOURCE_DIR})

#ADD_EXECUTABLE(gift SRC_LIST)
