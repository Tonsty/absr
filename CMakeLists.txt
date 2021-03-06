cmake_minimum_required(VERSION 2.8.12)

PROJECT(ABSR)
SET(ABSR_VERSION 1.1)

MESSAGE(STATUS "This is BINARY dir " ${PROJECT_BINARY_DIR})
MESSAGE(STATUS "This is SOURCE dir " ${PROJECT_SOURCE_DIR})

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
IF(CMAKE_HOST_UNIX)
	LIST(APPEND CMAKE_CXX_FLAGS -std=c++11)
ENDIF()
ADD_DEFINITIONS(-DABSR_PREINSTANTIATE)

IF(CMAKE_HOST_UNIX)
	# Set a default build type if none was specified
	IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	  MESSAGE(STATUS "Setting build type to 'Debug' as none was specified.")
	  SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
	  # Set the possible values of build type for cmake-gui
	  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
	    "MinSizeRel" "RelWithDebInfo")
	ENDIF()
ENDIF()

OPTION(WITH_QT4 "" ON)
OPTION(WITH_VTK "" ON)
OPTION(WITH_MATLAB "" ON)

# Find the Eigen library
FIND_PACKAGE(Eigen3 3.3.0 REQUIRED)

# Find the Qt4 library
IF(WITH_QT4)
	FIND_PACKAGE(Qt4 4.8.5 REQUIRED QtCore QtGui)
ENDIF(WITH_QT4)

# Find the VTK library
IF(WITH_VTK)
	FIND_PACKAGE(VTK 5.8 REQUIRED)
	include(${VTK_USE_FILE})
ENDIF(WITH_VTK)

IF(CMAKE_HOST_WIN32)
IF(WITH_MATLAB)
SET(MATLAB_INCLUDE_DIR "D:/Program Files/MATLAB/R2013a/extern/include")
SET(MATLAB_LIBARIES "D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libmat.lib"
					"D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libmx.lib"
					"D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libmex.lib"
					"D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libeng.lib")
ELSE(WITH_MATLAB)
SET(MATLAB_INCLUDE_DIR "")
SET(MATLAB_LIBARIES "")	
ENDIF(WITH_MATLAB)
ENDIF(CMAKE_HOST_WIN32)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
IF(CMAKE_HOST_WIN32)
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libstatic)
ELSEIF(CMAKE_HOST_UNIX)
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
ENDIF()

# uninstall target
CONFIGURE_FILE(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall
	COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR} 
	      			${MATLAB_INCLUDE_DIR}
	      			include)

SET(ABSR_INCLUDES include/tbsfitting.h
				  include/tbsfitting.hpp
				  include/boundingbox.h
				  include/typedefs.h
				  include/io.h
				  include/sdf.h
				  include/tensorbsplines.h
				  include/tensorbsplines.hpp
				  include/activetbs.h
				  include/activetbs.hpp
				  include/activetbsfitting.h
				  include/activetbsfitting.hpp
				  include/hierarchicaltbs.h
				  include/hierarchicaltbs.hpp
				  include/hierarchicaltbsfitting.h
				  include/hierarchicaltbsfitting.hpp
				  include/fastmarching.h
				  include/polygonizer.h
				  include/polygonizer.hpp
				  include/vtk_help.h
				  include/test.h
				  include/math_help.h) 
SET(ABSR_SOURCES src/main.cpp
				 src/boundingbox.cpp
				 src/io.cpp
				 src/sdf.cpp			 
				 src/fastmarching.cpp
				 src/polygonizer.c
				 src/polygonizer.cpp
				 src/tensorbsplines.cpp
				 src/test.cpp
				 src/vtk_help.cpp
				 src/math_help.cpp
				 src/tbsfitting.cpp
				 src/activetbs.cpp
				 src/activetbsfitting.cpp
				 src/hierarchicaltbs.cpp
				 src/hierarchicaltbsfitting.cpp)

ADD_EXECUTABLE(absr ${ABSR_INCLUDES} ${ABSR_SOURCES})
TARGET_LINK_LIBRARIES(absr ${VTK_LIBRARIES} ${MATLAB_LIBARIES})