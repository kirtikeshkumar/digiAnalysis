cmake_minimum_required(VERSION 2.8)
project(digiAnalysis)

# Added to detect compiler support for C++11
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

# Options for features
option(WAVES "Enable or Disable Waveform Analysis" OFF)
option(SMOOTH "Enable or Disable Waveform Smoothing" OFF)
option(CFD "Enable or Disable Waveform CFD" OFF)

if(WAVES)
    add_definitions(-DWAVES)
endif()
if(SMOOTH)
    add_definitions(-DSMOOTH)
endif()
if(CFD)
    add_definitions(-DCFD)
endif()

# Find ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

# Include directories
include_directories(inc)

# Gather source files from src directory
file(GLOB sources src/*.cpp)

# Gather test files from test directory
file(GLOB test_sources test/*.cpp)

# Loop through each test file and create an executable
foreach(test_file ${test_sources})
    # Get the filename without the directory path and extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create an executable for each test file
    add_executable(${test_name} ${test_file} ${sources})
    target_link_libraries(${test_name} ${ROOT_LIBRARIES})
endforeach()
