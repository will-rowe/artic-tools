######################################################################
# set up cmake
######################################################################
cmake_minimum_required(VERSION 3.12)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_DISABLE_SOURCE_CHANGES OFF)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# set a default build type if none was specified
set(COBS_DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${COBS_DEFAULT_BUILD_TYPE}' as none was specified")
  set(CMAKE_BUILD_TYPE "${COBS_DEFAULT_BUILD_TYPE}" CACHE
    STRING "Choose the type of build" FORCE)
  # set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

######################################################################
# build opts
######################################################################
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_TESTING "Build unit tests" ON)

######################################################################
# set up the project
######################################################################
project(ARTIC)

# set the variables
set(ARTIC_PROG_NAME artic-tools)
set(ARTIC_VERSION_MAJOR 0)
set(ARTIC_VERSION_MINOR 2)
set(ARTIC_VERSION_PATCH 2)
configure_file (
    "${PROJECT_SOURCE_DIR}/artic/version.hpp.in"
    "${PROJECT_SOURCE_DIR}/artic/version.hpp"
)

# create variables to collect some stuff during the build
set(ARTIC_INCLUDE_DIRS "")
set(ARTIC_LINK_LIBRARIES "")
set(ARTIC_DEFINITIONS "")

# print some messages
message(STATUS "Building from: CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")
message(STATUS "Installing to: CMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

######################################################################
# cpp settings
######################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CMAKE_COMPILER_IS_GNUCXX)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "7.0")
        message(FATAL_ERROR "GCC version must be at least 7.0!")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
        message(FATAL_ERROR "Clang version must be at least 9.0!")
    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wpedantic")

######################################################################
# check libraries
######################################################################
# pthreads
find_package(Threads REQUIRED)
set(ARTIC_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${ARTIC_LINK_LIBRARIES})
if(CMAKE_USE_PTHREADS_INIT)
  set(ARTIC_LINK_LIBRARIES pthread ${ARTIC_LINK_LIBRARIES})
endif()

# htslib
find_library(HTS_LIB hts HINTS /usr/local/lib/ /usr/lib/)
find_path(HTS_INCLUDE_DIR NAMES htslib/kseq.h HINTS "${HTSLIB_INSTALL_DIR}")
if (HTS_LIB AND HTS_INCLUDE_DIR)
    message(STATUS "Found htslib: ${HTS_LIB} and ${HTS_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "Could not find htslib")
endif()
set(ARTIC_INCLUDE_DIRS ${HTS_INCLUDE_DIR} ${ARTIC_INCLUDE_DIRS})
set(ARTIC_LINK_LIBRARIES ${HTS_LIB} ${ARTIC_LINK_LIBRARIES})

# boost (dynamic bitset)
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0 COMPONENTS filesystem REQUIRED) 

# cli11, spdlog
add_subdirectory(extlibs/CLI11)
add_subdirectory(extlibs/spdlog)

# rapidcsv
include_directories(${PROJECT_SOURCE_DIR}/extlibs/rapidcsv/src)

# google test
if(BUILD_TESTING)
    add_subdirectory(extlibs/googletest)
    include(GoogleTest)
endif(BUILD_TESTING)

######################################################################
# build
######################################################################
add_subdirectory(artic)
add_subdirectory(app)

######################################################################
# test
######################################################################
if(BUILD_TESTING)
    message(STATUS "Running tests")
    enable_testing()
    add_subdirectory(tests)
endif(BUILD_TESTING)