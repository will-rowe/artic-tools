file(COPY data DESTINATION .)

file(
  GLOB_RECURSE TEST_FILES
  "${PROJECT_SOURCE_DIR}/tests/*.[ch]pp"
)

add_executable(${PROJECT_NAME}_tests ${TEST_FILES})
target_link_libraries(
  ${PROJECT_NAME}_tests 
  artic_static 
  gtest 
  gtest_main
)
target_compile_definitions(
  ${PROJECT_NAME}_tests 
  PRIVATE 
  TEST_DATA_PATH="${PROJECT_SOURCE_DIR}/tests/data/"
)
gtest_add_tests(${PROJECT_NAME}_tests "" ${TEST_FILES})
