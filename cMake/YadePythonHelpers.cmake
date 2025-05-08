##########################################################################
#  Copyright (C) 2019 by Kneib Francois                                  #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################


# Inspired by http://www.cmake.org/pipermail/cmake/2011-January/041666.html
#
# - Find Python Module
FUNCTION(FIND_PYTHON_MODULE module)
  #Set a variable with the module name in upper case
  STRING(TOUPPER ${module} module_upper)

  #Reset result value as this function can be called multiple times while testing different python versions:
  UNSET(PY_${module_upper} CACHE)
  UNSET(${module_upper}_INCLUDE_DIR CACHE)

  IF(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
    SET(${module_upper}_FIND_REQUIRED TRUE)
  ENDIF(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")

  EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    #Since tkinter does not have __version__ attribute, use TkVersion attribute instead to get the version
    "import re, ${module}; \
    location = re.compile('/__init__.py.*').sub('', ${module}.__file__); \
    include_dir = ${module}.get_include() if hasattr(${module}, 'get_include') else None; \
    version = ${module}.TkVersion if hasattr(${module}, 'TkVersion') else ${module}.__version__; \
    print(location, include_dir, version);"
    RESULT_VARIABLE _${module}_status
    ERROR_VARIABLE _${module}_error
    OUTPUT_VARIABLE _${module}_output
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

  IF(_${module}_status MATCHES 0)
    #Split the _${module}_output into a list
    STRING(REPLACE " " ";" _${module}_output_list ${_${module}_output})

    #Get location from the first element of the list
    LIST(GET _${module}_output_list 0 _${module}_location)

    #Get include dir from the second element of the list
    LIST(GET _${module}_output_list 1 _${module}_include_dir)

    #Set MODULE_FOUND variable
    IF(_${module}_include_dir MATCHES "None")
      SET(PY_${module_upper} ${_${module}_location} CACHE STRING "Location of Python module ${module}")
      FIND_PACKAGE_HANDLE_STANDARD_ARGS(${module_upper} DEFAULT_MSG PY_${module_upper})
    ELSEIF(_${module}_include_dir MATCHES "Traceback")  # Safeguard for mpi4py
      SET(${module_upper}_FOUND False)
      IF(${module_upper}_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "${module} include directory not found : matches Traceback")
        RETURN()
      ELSE()
        MESSAGE("${module} include directory not found : matches Traceback")
      ENDIF()
    ELSE()  # for numpy, mpi4py or other libraries with get_include() method
      SET(${module_upper}_INCLUDE_DIR ${_${module}_include_dir} CACHE STRING "${module} include directory")
      FIND_PACKAGE_HANDLE_STANDARD_ARGS(${module_upper} DEFAULT_MSG ${module_upper}_INCLUDE_DIR)
      INCLUDE_DIRECTORIES(${${module_upper}_INCLUDE_DIR})
    ENDIF()

    #Concatenate the rest of the list to create the version in case there is space in the python __version__
    #Can't use LIST(SUBLIST ...) as it is not compatible with old versions of cmake
    SET(_${module}_version)
    LIST(LENGTH _${module}_output_list len)
    MATH(EXPR end_index "${len} - 1")
    FOREACH(i RANGE 2 ${end_index})
      LIST(GET _${module}_output_list ${i} item)
      LIST(APPEND _${module}_version ${item})
    ENDFOREACH()
    STRING(REPLACE ";" " " ${module_upper}_VERSION ${_${module}_version})

    #Get the version major, minor and patch depending of the version format
    #Feel free to add other formats, for now, only x.x.x, x.x or (x, x) formats are supported
    IF(${module_upper}_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)$")
      SET(${module_upper}_VERSION_MAJOR "${CMAKE_MATCH_1}" )
      SET(${module_upper}_VERSION_MINOR "${CMAKE_MATCH_2}")
      SET(${module_upper}_VERSION_PATCH "${CMAKE_MATCH_3}")
    ELSEIF(${module_upper}_VERSION MATCHES "^([0-9]+)\\.([0-9]+)$")
      SET(${module_upper}_VERSION_MAJOR "${CMAKE_MATCH_1}")
      SET(${module_upper}_VERSION_MINOR "${CMAKE_MATCH_2}")
      SET(${module_upper}_VERSION_PATCH "0")
    ELSEIF(${module_upper}_VERSION MATCHES "^\\(([0-9]+),([0-9]+)\\)$")
      SET(${module_upper}_VERSION_MAJOR "${CMAKE_MATCH_1}")
      SET(${module_upper}_VERSION_MINOR "${CMAKE_MATCH_2}")
      SET(${module_upper}_VERSION_PATCH "0")
    ENDIF()

  ELSE(_${module}_status MATCHES 0)
    SET(${module_upper}_FOUND FALSE)
    IF(${module_upper}_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "${module} import failure:\n${_${module}_error}")
      RETURN()
    ENDIF()
  ENDIF(_${module}_status MATCHES 0)

  #Set the variables with PARENT_SCOPE in order to access them outside of the function
  SET(${module_upper}_FOUND ${${module_upper}_FOUND} PARENT_SCOPE)
  IF(${module_upper}_FOUND)
    SET(${module_upper}_VERSION ${${module_upper}_VERSION} PARENT_SCOPE)
    SET(${module_upper}_VERSION_MAJOR ${${module_upper}_VERSION_MAJOR} PARENT_SCOPE)
    SET(${module_upper}_VERSION_MINOR ${${module_upper}_VERSION_MINOR} PARENT_SCOPE)
    SET(${module_upper}_VERSION_PATCH ${${module_upper}_VERSION_PATCH} PARENT_SCOPE)
  ELSE()
    UNSET(${module_upper}_VERSION PARENT_SCOPE)
    UNSET(${module_upper}_VERSION_MAJOR PARENT_SCOPE)
    UNSET(${module_upper}_VERSION_MINOR PARENT_SCOPE)
    UNSET(${module_upper}_VERSION_PATCH PARENT_SCOPE)
  ENDIF()
ENDFUNCTION(FIND_PYTHON_MODULE)


# Find PythonLibs, all Python packages needed by yade, and libboost-python. Must be used after a call to FIND_PACKAGE(PythonInterp)
FUNCTION(FIND_PYTHON_PACKAGES)
	SET(ALL_PYTHON_DEPENDENCIES_FOUND FALSE PARENT_SCOPE)
	SET(fail_message "Failed to import dependencies for Python version ${PYTHON_VERSION_STRING}. NOT FOUND:")

	UNSET(PYTHON_LIBRARY CACHE)
	UNSET(PYTHON_INCLUDE_DIR CACHE)
	FIND_PACKAGE(PythonLibs QUIET)
	IF(NOT PYTHONLIBS_FOUND)
		MESSAGE(${fail_message} PythonLibs)
		RETURN()
	ENDIF()

	# BEGIN find Boost for py_version
	IF ( NOT LocalBoost )
		SET(LocalBoost "1.47.0") # Minimal required Boost version
	ENDIF ( NOT LocalBoost )
	# Next loop is due to libboost-pythonXXX naming mismatch between ubuntu versions and debian versions, so try three possibilities that cover all distros.
	FOREACH(PYTHON_PREFIX python python-py python${PYTHON_VERSION_MAJOR}-py) #boost>1.67 should pick-up the first one.
		IF(ENABLE_LOGGER)
			FIND_PACKAGE(Boost ${LocalBoost}  QUIET COMPONENTS ${PYTHON_PREFIX}${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} thread filesystem iostreams regex serialization system date_time log)
		ELSE(ENABLE_LOGGER)
			FIND_PACKAGE(Boost ${LocalBoost}  QUIET COMPONENTS ${PYTHON_PREFIX}${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} thread filesystem iostreams regex serialization system date_time)
		ENDIF(ENABLE_LOGGER)
		IF(Boost_FOUND)
			# for some reason boost_python37 is found but not linked with boost 1.71, we add it here (is it a specific issue within NIX?)
			IF (${Boost_VERSION} GREATER 107100 OR ${Boost_VERSION} EQUAL 107100) #maybe it should start at boost 1.67?
				MESSAGE("Boost_VERSION=${Boost_VERSION}, adding boost_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} lib")
				SET(Boost_LIBRARIES "${Boost_LIBRARIES};libboost_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.so")

			ENDIF()
			BREAK()
		ENDIF()
	ENDFOREACH()

	IF(NOT Boost_FOUND) # for opensuze
		IF(ENABLE_LOGGER)
			FIND_PACKAGE(Boost ${LocalBoost}  QUIET COMPONENTS python-py${PYTHON_VERSION_MAJOR} thread filesystem iostreams regex serialization system date_time log)
		ELSE(ENABLE_LOGGER)
			FIND_PACKAGE(Boost ${LocalBoost}  QUIET COMPONENTS python-py${PYTHON_VERSION_MAJOR} thread filesystem iostreams regex serialization system date_time)
		ENDIF(ENABLE_LOGGER)
	ENDIF()

	IF(NOT Boost_FOUND) #as we try multiple python prefixes we have to handle manually the required behavior: fail if we didn't found boost
		MESSAGE(${fail_message} libboost-python)
		RETURN()
	ENDIF()
	# END find Boost for py_version

	# Find Python modules and set the version variable in the parent scope which is CMakeLists.txt
	FOREACH(module IN ITEMS IPython numpy matplotlib pygraphviz Xlib sphinx tkinter)
		IF(${module} MATCHES "tkinter" AND ${PYTHON_VERSION_MAJOR} EQUAL 2)
			SET(module "Tkinter")
		ENDIF()
		FIND_PYTHON_MODULE(${module} REQUIRED)
		STRING(TOUPPER ${module} module_upper)
		SET(${module_upper}_FOUND ${${module_upper}_FOUND} PARENT_SCOPE)
		IF(${module_upper}_FOUND)
			MESSAGE(STATUS "${module} version found: ${${module_upper}_VERSION}")
			SET(${module_upper}_VERSION ${${module_upper}_VERSION} PARENT_SCOPE)
			SET(${module_upper}_VERSION_MAJOR ${${module_upper}_VERSION_MAJOR} PARENT_SCOPE)
			SET(${module_upper}_VERSION_MINOR ${${module_upper}_VERSION_MINOR} PARENT_SCOPE)
			SET(${module_upper}_VERSION_PATCH ${${module_upper}_VERSION_PATCH} PARENT_SCOPE)
		ELSE()
			MESSAGE(${fail_message} ${module})
			UNSET(${module_upper}_VERSION PARENT_SCOPE)
			UNSET(${module_upper}_VERSION_MAJOR PARENT_SCOPE)
			UNSET(${module_upper}_VERSION_MINOR PARENT_SCOPE)
			UNSET(${module_upper}_VERSION_PATCH PARENT_SCOPE)
			RETURN()
		ENDIF()
	ENDFOREACH()

	# NOTE: If we are here, we found a suitable Python version with all packages needed.
	SET(ALL_PYTHON_DEPENDENCIES_FOUND TRUE PARENT_SCOPE)
	#Export findpythonlibs vars to global parent scope:
	FOREACH(pythonlibs_var PYTHONLIBS_FOUND PYTHON_LIBRARIES PYTHON_INCLUDE_PATH PYTHON_INCLUDE_DIRS PYTHONLIBS_VERSION_STRING NUMPY_VERSION_MAJOR NUMPY_VERSION_MINOR)
		SET(${pythonlibs_var} ${${pythonlibs_var}} PARENT_SCOPE)
	ENDFOREACH()
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
	#Export findboost vars to global parent scope:
	FOREACH(boost_var boost_FOUND Boost_INCLUDE_DIRS Boost_LIBRARY_DIRS Boost_LIBRARIES Boost_<C>_FOUND Boost_<C>_LIBRARY Boost_VERSION Boost_LIB_VERSION Boost_MAJOR_VERSION Boost_MINOR_VERSION Boost_SUBMINOR_VERSION)
		SET(${boost_var} ${${boost_var}} PARENT_SCOPE)
	ENDFOREACH()
	# for checking purpose
	MESSAGE("--   Boost_VERSION: " ${Boost_VERSION})
	MESSAGE("--   Boost_LIB_VERSION: " ${Boost_LIB_VERSION})
	MESSAGE("--   Boost_INCLUDE_DIRS: " ${Boost_INCLUDE_DIRS})
	MESSAGE("--   Boost_LIBRARIES: " ${Boost_LIBRARIES})

ENDFUNCTION(FIND_PYTHON_PACKAGES)

# Did findpythoninterp found the python version we want ? Output in PYTHON_VERSION_MATCH.
FUNCTION(PYTHON_VERSION_MATCHES version_number)
	SET(PYTHON_VERSION_MATCH FALSE PARENT_SCOPE)
	string(REGEX MATCH "([0-9]+)\\.([0-9]+)" _ ${version_number})
	set(ver_major ${CMAKE_MATCH_1})
	set(ver_minor ${CMAKE_MATCH_2})
	MESSAGE("Trying python version: " ${version_number} " parsed as " ${ver_major} " " ${ver_minor})

	IF(NOT (${PYTHON_VERSION_MAJOR} EQUAL ${ver_major}))
		RETURN()
	ENDIF()

	IF(NOT(${PYTHON_VERSION_MINOR} EQUAL ${ver_minor}))
		RETURN()
	ENDIF()
	#if we are here we match major and minor
	SET(PYTHON_VERSION_MATCH TRUE PARENT_SCOPE)
ENDFUNCTION(PYTHON_VERSION_MATCHES version_number)
