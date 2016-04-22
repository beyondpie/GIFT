FIND_PACKAGE(Boost COMPONENTS program_options thread system exception)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "The boost lib file dir is " ${Boost_LIBRRIES})

  MESSAGE(STATUS "gift binary dir is " ${PROJECT_BINARY_DIR})
  MESSAGE(STATUS "gift source dir is " ${PROJECT_SOURCE_DIR})
  ADD_EXECUTABLE(gift
    gift.hpp
#    globalvar_gift.hpp
#    EMGift.hpp
    EMGift.cpp
#    parameterGift.hpp
    parameterGift.cpp
    functionsGift.cpp
    main.cpp)
  TARGET_LINK_LIBRARIES(gift ${Boost_LIBRARIES})
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Uable to find correct Boost library. Did you set BOOST_ROOT?")
ENDIF(Boost_FOUND)
