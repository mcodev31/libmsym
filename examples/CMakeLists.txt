cmake_minimum_required (VERSION 2.8.11)

add_executable (msym_example msym_example.c example.c)
add_executable (msym_tex tex_character_table)

target_link_libraries (msym_example LINK_PUBLIC msym)
target_link_libraries (msym_tex LINK_PUBLIC msym)
