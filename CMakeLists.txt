cmake_minimum_required(VERSION 2.8.12)

PROJECT(ABSR)
SET(ABSR_VERSION 1.1)

MESSAGE(STATUS "This is BINARY dir " ${PROJECT_BINARY_DIR})
MESSAGE(STATUS "This is SOURCE dir " ${PROJECT_SOURCE_DIR})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

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

IF(WITH_MATLAB)
SET(MATLAB_INCLUDE_DIR "D:/Program Files/MATLAB/R2013a/extern/include")
SET(MATLAB_LIBARIES "D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libmat.lib"
					"D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libmx.lib"
					"D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libmex.lib"
					"D:/Program Files/MATLAB/R2013a/extern/lib/win64/microsoft/libeng.lib")
ENDIF(WITH_MATLAB)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
IF(CMAKE_HOST_WIN32)
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libstatic)
ELSEIF(CMAKE_HOST_WIN32)
	SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
ENDIF(CMAKE_HOST_WIN32)

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

SET(ABSR_INCLUDES include/absr.h
				  include/boundingbox.h
				  include/bspline.h
				  include/common.h
				  include/io.h
				  include/sdf
				  include/tensorbsplines.h) 
SET(ABSR_SOURCES src/absr.cpp 
				 src/main.cpp
				 src/boundingbox.cpp
				 src/bspline.cpp
				 src/io.cpp
				 src/sdf.cpp
				 src/tensorbsplines.cpp)

ADD_EXECUTABLE(absr ${ABSR_INCLUDES} ${ABSR_SOURCES})
TARGET_LINK_LIBRARIES(absr ${VTK_LIBRARIES} ${MATLAB_LIBARIES})