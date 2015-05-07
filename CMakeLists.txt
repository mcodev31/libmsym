cmake_minimum_required (VERSION 2.8.11)

file (GLOB source src/*.c)
add_library (msym ${source})
if (BUILD_SHARED_LIBS)
        set_target_properties (msym PROPERTIES
                POSITION_INDEPENDENT_CODE "TRUE")
endif ()
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
