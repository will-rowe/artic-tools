######################################################################
# set up cmake
######################################################################
cmake_minimum_required(VERSION 3.10)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_DISABLE_SOURCE_CHANGES OFF)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# set a default build type if none was specified
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
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
project(ARTIC LANGUAGES CXX)

# set the variables
set(ARTIC_PROG_NAME artic-tools)
set(ARTIC_VERSION_MAJOR 0)
set(ARTIC_VERSION_MINOR 3)
set(ARTIC_VERSION_PATCH 0)
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
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wpedantic -O0 -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -O3")
endif()

######################################################################
# check libraries
######################################################################
# pthreads
find_package(Threads REQUIRED)

# boost (dynamic bitset, filesystem, json)
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0 COMPONENTS filesystem REQUIRED)

# curl
set(CURL_LIBRARY "-lcurl") 
find_package(CURL REQUIRED) 

# zlib
find_package(ZLIB REQUIRED)

# link libraries
set(ARTIC_LINK_LIBRARIES ZLIB::ZLIB ${Boost_LIBRARIES} ${CURL_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${ARTIC_LINK_LIBRARIES})
if(CMAKE_USE_PTHREADS_INIT)
  set(ARTIC_LINK_LIBRARIES pthread ${ARTIC_LINK_LIBRARIES})
endif()

# external dependencies
include(${PROJECT_SOURCE_DIR}/extlibs/external.cmake)

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
