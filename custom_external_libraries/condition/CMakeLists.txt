set(CONDITION_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/r8lib.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/condition.cpp
)

add_library(condition STATIC ${CONDITION_SOURCE})
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(condition PUBLIC /wd4804)
endif()
