cmake_minimum_required (VERSION 2.8.11)

#project(libmsym C) # unfortunately GenerateExportHeader doesn't support C projects, which sort of makes this project reqire c++ =/

project(libmsym)

option(MSYM_BUILD_EXAMPLES "Build example executables" OFF)
option(MSYM_BUILD_PYTHON "Build python binding" OFF)

include (GenerateExportHeader)

set(LIBMSYM_MAJOR_VERSION 0)
set(LIBMSYM_MINOR_VERSION 2)
set(LIBMSYM_PATCH_VERSION 4)
set(LIBMSYM_VERSION ${LIBMSYM_MAJOR_VERSION}.${LIBMSYM_MINOR_VERSION}.${LIBMSYM_PATCH_VERSION})

set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for msym library")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/libmsym)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

include_directories("${PROJECT_SOURCE_DIR}")

file (GLOB source src/*.c)
add_library (msym ${source})

if (BUILD_SHARED_LIBS)
	set_target_properties(msym PROPERTIES
		POSITION_INDEPENDENT_CODE TRUE
		C_VISIBILITY_PRESET hidden
		VISIBILITY_INLINES_HIDDEN 1)

else()
	set_target_properties(msym PROPERTIES COMPILE_FLAGS -DMSYM_EXPORTS_BUILT_AS_STATIC)
endif ()

include_directories("${PROJECT_BINARY_DIR}")

set_target_properties(msym PROPERTIES VERSION ${LIBMSYM_VERSION})
set_target_properties(msym PROPERTIES SOVERSION ${LIBMSYM_MAJOR_VERSION}.${LIBMSYM_MINOR_VERSION})

generate_export_header(msym
             BASE_NAME msym
             EXPORT_MACRO_NAME MSYM_EXPORT
             EXPORT_FILE_NAME msym_EXPORTS.h
	     STATIC_DEFINE MSYM_EXPORTS_BUILT_AS_STATIC)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
#set_target_properties (msym PROPERTIES
#        C_STANDARD 99
#        C_STANDARD_REQUIRED ON)

#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -flto")
#set (CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -flto" )
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fsanitize=address -fno-omit-frame-pointer")
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fsanitize=memory -fno-omit-frame-pointer -fsanitize-memory-track-origins")

set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -DLIBMSYM_DEBUG")

set_target_properties(msym PROPERTIES PUBLIC_HEADER "src/msym.h;src/msym_error.h")

if(UNIX)
        target_link_libraries(msym m)
endif(UNIX)

export(TARGETS msym FILE "${PROJECT_BINARY_DIR}/libmsymTargets.cmake")

export(PACKAGE libmsym)

file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(libmsymConfig.cmake.in "${PROJECT_BINARY_DIR}/libmsymConfig.cmake" @ONLY)
set(CONF_INCLUDE_DIRS "\${LIBMSYM_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(libmsymConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libmsymConfig.cmake" @ONLY)
configure_file(libmsymConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/libmsymConfigVersion.cmake" @ONLY)

install(FILES
    "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libmsymConfig.cmake"
    "${PROJECT_BINARY_DIR}/libmsymConfigVersion.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}")

# EXPORT HEADER
install(FILES "${PROJECT_BINARY_DIR}/msym_EXPORTS.h" DESTINATION "${INSTALL_INCLUDE_DIR}/libmsym")

install(TARGETS msym
        EXPORT libmsymTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
        ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
        PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/libmsym"
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}" )

install(EXPORT libmsymTargets DESTINATION "${INSTALL_CMAKE_DIR}")

if(MSYM_BUILD_EXAMPLES)
	add_subdirectory(examples)
endif()

if(MSYM_BUILD_PYTHON)
	if(BUILD_SHARED_LIBS)
		add_subdirectory(bindings/python)
	else()
		message(SEND_ERROR "python requires libmsym to be built as a shared library, use -DBUILD_SHARED_LIBS:BOOL=ON")
	endif()
endif()
