cmake_minimum_required (VERSION 2.8.11)

find_program(PYTHON "python")
if (PYTHON)
	set(MSYM_PYTHON_PACKAGE_DIR  "${CMAKE_CURRENT_BINARY_DIR}/libmsym")
	set(MSYM_PYTHON_SETUP_IN     "${CMAKE_CURRENT_SOURCE_DIR}/cmake/setup.py.in")
	set(MSYM_PYTHON_INIT_C_IN    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/__init__.py.in")
	set(MSYM_PYTHON_CODE_IN      "${CMAKE_CURRENT_SOURCE_DIR}/libmsym/libmsym.py")
	set(MSYM_PYTHON_SETUP        "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
	set(MSYM_PYTHON_INIT_G_IN    "${CMAKE_CURRENT_BINARY_DIR}/__init__.py.in")
	set(MSYM_PYTHON_INIT         "${MSYM_PYTHON_PACKAGE_DIR}/__init__.py")
	set(MSYM_PYTHON_CODE         "${MSYM_PYTHON_PACKAGE_DIR}/libmsym.py")
	set(MSYM_PYTHON_STAMP        "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")
	set(MSYM_PYTHON_INSTALL_OPTS "" CACHE STRING "install options for setup.py")

	# There simply has to be a better way to get distutils to install a shared library
	# from a custom path without overriding all of the command classes,
	# get cmake to tell me where the library is located in a less convoluted way
	# or get python to use other means than ldconfig to load shared libraries... this is ridiculous
	configure_file(${MSYM_PYTHON_SETUP_IN} ${MSYM_PYTHON_SETUP})
	configure_file(${MSYM_PYTHON_INIT_C_IN} ${MSYM_PYTHON_INIT_G_IN})
	file(GENERATE OUTPUT ${MSYM_PYTHON_INIT} INPUT ${MSYM_PYTHON_INIT_G_IN})
	file(COPY ${MSYM_PYTHON_CODE_IN} DESTINATION ${MSYM_PYTHON_PACKAGE_DIR}) 

	add_custom_command(OUTPUT ${MSYM_PYTHON_STAMP}
			COMMAND ${PYTHON} ${MSYM_PYTHON_SETUP} build
			COMMAND ${CMAKE_COMMAND} -E touch ${MSYM_PYTHON_STAMP}
			DEPENDS ${MSYM_PYTHON_INIT_IN} ${MSYM_PYTHON_SETUP_IN} ${MSYM_PYTHON_CODE_IN})

	add_custom_target(target ALL DEPENDS ${MSYM_PYTHON_STAMP})

	install(CODE "execute_process(COMMAND ${PYTHON} ${MSYM_PYTHON_SETUP} install ${MSYM_PYTHON_INSTALL_OPTS})")
endif()

