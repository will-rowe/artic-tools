file(COPY data DESTINATION .)

file(
  GLOB_RECURSE TEST_FILES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.[ch]pp
)

add_executable(artic_tests ${TEST_FILES})
target_compile_definitions(artic_tests PRIVATE artic_test)
target_link_libraries(artic_tests artic_static gtest gtest_main)
gtest_add_tests(artic_tests "" ${TEST_FILES})
