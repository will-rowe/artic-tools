# libs
## htslib
find_library(HTS_LIB hts HINTS /usr/local/lib/ /usr/lib/)
if (HTS_LIB)
    message(STATUS "found htslib: ${HTS_LIB}")
else()
    message(STATUS "could not find htslib, attempting to compile")
    set(HTSLIB_SOURCE_DIR ${PROJECT_SOURCE_DIR}/extlibs/htslib)
    execute_process(
        COMMAND git submodule update --init --recursive -- ${HTSLIB_SOURCE_DIR}
        WORKING_DIRECTORY ${HTSLIB_SOURCE_DIR}
        COMMAND autoheader && autoconf && sudo make && sudo make install
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )    
endif()

## cli11
if(NOT TARGET CLI11::CLI11)
  set(CLI11_SOURCE_DIR ${PROJECT_SOURCE_DIR}/extlibs/CLI11)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${CLI11_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_subdirectory(${CLI11_SOURCE_DIR})
endif()

## kseq++
if(NOT TARGET kseq++::kseq++)
  set(KSEQPP_SOURCE_DIR ${PROJECT_SOURCE_DIR}/extlibs/kseq++)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${KSEQPP_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_subdirectory(${KSEQPP_SOURCE_DIR})
endif()

## spdlog
if(NOT TARGET spdlog::spdlog)
  set(SPDLOG_SOURCE_DIR ${PROJECT_SOURCE_DIR}/extlibs/spdlog)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${SPDLOG_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  add_subdirectory(${SPDLOG_SOURCE_DIR})
endif()

## google test
if(BUILD_TESTING)
    if(NOT TARGET gtest::gtest)
        set(GTEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/extlibs/googletest)
        execute_process(
            COMMAND git submodule update --init --recursive -- ${GTEST_SOURCE_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_subdirectory(${PROJECT_SOURCE_DIR}/extlibs/googletest)
        include(GoogleTest)
    endif()
endif(BUILD_TESTING)

## link
set(ARTIC_LINK_LIBRARIES ${HTS_LIB} CLI11::CLI11 kseq++::kseq++ spdlog::spdlog ${ARTIC_LINK_LIBRARIES})

# includes
## htslib
find_path(HTS_INCLUDE_DIR NAMES htslib/kseq.h HINTS "${HTS_LIB}")
if (NOT HTS_INCLUDE_DIR)
    message(FATAL_ERROR "could not include htslib")
endif()

## rapidcsv
set(RAPIDCSV ${PROJECT_SOURCE_DIR}/extlibs/rapidcsv)
if (NOT EXISTS ${RAPIDCSV})
    execute_process(
        COMMAND git submodule update --init --recursive -- ${RAPIDCSV}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
set(RAPIDCSV_INCLUDE_DIR ${RAPIDCSV}/src)

## flat_hash_map
set(FLATHASHMAP ${PROJECT_SOURCE_DIR}/extlibs/flat_hash_map)
if (NOT EXISTS ${FLATHASHMAP})
    execute_process(
        COMMAND git submodule update --init --recursive -- ${FLATHASHMAP}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

## add them
set(ARTIC_INCLUDE_DIRS ${HTS_INCLUDE_DIR} ${RAPIDCSV_INCLUDE_DIR} ${FLATHASHMAP} ${ARTIC_INCLUDE_DIRS})
