cmake_minimum_required (VERSION 2.8.11)

include_directories(${PROJECT_BINARY_DIR} "${PROJECT_SOURCE_DIR}/src")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
#set_target_properties (msym PROPERTIES
#        C_STANDARD 99
#        C_STANDARD_REQUIRED ON)


#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -flto")
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fsanitize=address -fno-omit-frame-pointer")
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fsanitize=memory -fno-omit-frame-pointer -fsanitize-memory-track-origins")
#set (CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -flto" )

add_executable (msym_example msym_example.c example.c)
add_executable (msym_tex tex_character_table)

target_link_libraries (msym_example LINK_PUBLIC msym)
target_link_libraries (msym_tex LINK_PUBLIC msym)
