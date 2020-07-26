file(COPY data DESTINATION .)

file(
  GLOB_RECURSE TEST_FILES
  "${PROJECT_SOURCE_DIR}/tests/*.[ch]pp"
)

add_executable(artic_tests ${TEST_FILES})

target_link_libraries(
  artic_tests 
  artic_static 
  gtest 
  gtest_main
)
target_compile_definitions(
  artic_tests 
  PRIVATE 
  TEST_DATA_PATH="${PROJECT_SOURCE_DIR}/tests/data/"
)

gtest_add_tests(artic_tests "" ${TEST_FILES})
