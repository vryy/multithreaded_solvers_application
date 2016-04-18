set(CMAKE_INCLUDE_CURRENT_DIR ON)

#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/SRC )

set( HSL_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/ddeps.f
        ${CMAKE_CURRENT_SOURCE_DIR}/mc29d.f
        ${CMAKE_CURRENT_SOURCE_DIR}/mc75d.f
        ${CMAKE_CURRENT_SOURCE_DIR}/hsl.cpp
     )

###############################################################
add_library(hsl STATIC ${HSL_SRC})
# install(TARGETS hsl DESTINATION libs )
