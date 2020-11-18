file(
    GLOB ARTIC_SOURCES
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/*.[ch]pp
)

add_library(
    artic_static
    STATIC
    ${ARTIC_SOURCES}
)

target_compile_definitions(artic_static PUBLIC ${ARTIC_DEFINITIONS})
target_include_directories(artic_static PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(artic_static SYSTEM PUBLIC ${ARTIC_INCLUDE_DIRS})
target_link_libraries(artic_static PUBLIC ${ARTIC_LINK_LIBRARIES})
