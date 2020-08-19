file(
  GLOB_RECURSE EXECUTABLE_FILES
  ${PROJECT_SOURCE_DIR}/app/*.cpp
)

foreach(target ${EXECUTABLE_FILES})
  get_filename_component(target_name ${target} NAME_WE)
  add_executable(${target_name} ${target})
  target_link_libraries(${target_name} artic_static CLI11::CLI11 ${Boost_LIBRARIES})
endforeach(target)