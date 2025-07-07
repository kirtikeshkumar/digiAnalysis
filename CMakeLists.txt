cmake_minimum_required(VERSION 2.8)
project(MST)

#Added to detect compiler support for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
elseif(COMPILER_SUPPORTS_CXX0X)
message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()


option(WAVES "Enable or Disable Waveform Analysis" ON)
option(SMOOTH "Enable or Disable Waveform Smoothing" OFF)
option(CFD "Enable or Disable Waveform CFD" OFF)

if(WAVES)
add_definitions(-DWAVES)
endif(WAVES)

if(SMOOTH)
add_definitions(-DSMOOTH)
endif(SMOOTH)

if(CFD)
add_definitions(-DCFD)
endif(CFD)

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
#message(STATUS ${ROOT_LIBRARIES})

include_directories(inc)

file(GLOB sources src/*.cpp)

#add_executable(Test_Dummy test/Test_Dummy.cpp ${sources})
#target_link_libraries(Test_Dummy ${ROOT_LIBRARIES})

SET(TESTSDIR ${TESTSDIR} ${CMAKE_SOURCE_DIR}/test)

#COLLECTING ALL TEST FILES
macro(collect_tests TESTS)
  foreach(TEST ${TESTS})
    file(GLOB tests ${tests} ${TEST}/*.cpp)
  endforeach()
endmacro()

macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE} ${headers} ${sources})
    target_link_libraries(${TARGET_NAME} ${ROOT_LIBRARIES})
 endforeach()
endmacro()

collect_tests("${TESTSDIR}")
build_executables("${tests}")

