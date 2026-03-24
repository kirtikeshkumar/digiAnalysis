cmake_minimum_required(VERSION 3.5)
project(digiAnalysis)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
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
option(WAVES "Enable or Disable Waveform Analysis" ON)
option(SMOOTH "Enable or Disable Waveform Smoothing" OFF)
option(CFD "Enable or Disable Waveform CFD" OFF)

if(WAVES)
    add_compile_definitions(WAVES)
endif()
if(SMOOTH)
    add_compile_definitions(SMOOTH)
endif()
if(CFD)
    add_compile_definitions(CFD)
endif()

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Spectrum)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
message(STATUS "ROOT_USE_FILE = ${ROOT_USE_FILE}")
message(STATUS "ROOT_INCLUDE_DIRS = ${ROOT_INCLUDE_DIRS}")

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)

set(DICT_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/WaveForm.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/singleHits.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/Pair.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/PSBar.h
)


ROOT_GENERATE_DICTIONARY(G__digiAnalysis
    ${DICT_HEADERS}
    LINKDEF inc/DictLinkDef.h 
    MODULE digiAnalysisLib
)


# Gather source files from src directory
file(GLOB sources src/*.cpp)

add_library(digiAnalysisLib SHARED
    ${sources}
    G__digiAnalysis.cxx
)

target_link_libraries(digiAnalysisLib ${ROOT_LIBRARIES} ROOT::Spectrum)
target_include_directories(digiAnalysisLib PUBLIC inc)

# Gather test files from test directory
file(GLOB test_sources CONFIGURE_DEPENDS test/*.cpp)

# Loop through each test file and create an executable
foreach(test_file ${test_sources})
    # Get the filename without the directory path and extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create an executable for each test file
    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} digiAnalysisLib ${ROOT_LIBRARIES} ROOT::Spectrum)
    #message(STATUS "Test sources: ${test_name}")
endforeach()
