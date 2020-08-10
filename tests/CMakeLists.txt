file(
  GLOB_RECURSE TEST_FILES
  "${PROJECT_SOURCE_DIR}/tests/*.[ch]pp"
)

add_executable(artic-tools_tests ${TEST_FILES})
  
target_link_libraries(
  artic-tools_tests
  artic_static 
  gtest
  gtest_main
)
    
file(COPY data DESTINATION .)
target_compile_definitions(
  artic-tools_tests 
  PRIVATE 
  TEST_DATA_PATH="${PROJECT_SOURCE_DIR}/tests/data/"
)

gtest_add_tests(TARGET artic-tools_tests)
