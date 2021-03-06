# MIT License

# Copyright (c) 2018 Benjamin Bercovici

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



# @file   CMakeLists.txt
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   July, 2017
# @brief  CMake listing enabling compilation and installation of the CGAL Interface library


################################################################################
#
# 								User-defined paths
#						Should be checked for consistency
#						Before running 'cmake ..' in build dir
#
################################################################################

# OMP-friendly GCC Major version (will only be used on Mac OS systems)
set (OMP_FRIENDLY_GCC_MAJOR_VERSION 7)

# Path to OMP-Friendly GCC Compiler (will only be used on Mac OS systems)
# set (OMP_FRIENDLY_GCC_PATH /usr/local/Cellar/gcc/7.3.0_1/bin/)

################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################

cmake_minimum_required(VERSION 3.0.0)


if (EXISTS /home/bebe0705/.am_fortuna)
	set(IS_FORTUNA ON)
	message("-- This is Fortuna")

else()
	set(IS_FORTUNA OFF)
	message("-- This is not Fortuna")
endif()

if (${USE_GCC})
	include(cmake/FindOmpGcc.cmake)
else()
	set(CMAKE_C_COMPILER /usr/bin/gcc CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER /usr/bin/g++ CACHE STRING "C++ Compiler" FORCE)
endif()


# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME CGAL_interface CACHE STRING "Name of library to be created.")

project(${LIB_NAME})


# Specify the version used
if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()


cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)

# Enable C++17 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD 11)

# Compiler flags
add_definitions(-Wall -O2)

# Include directories
include_directories(include)

# Find Boost
find_package(Boost COMPONENTS filesystem system REQUIRED) 
include_directories(${Boost_INCLUDE_DIRS}) 


# Find CGAL
find_package(CGAL REQUIRED)
include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

# Find Eigen
find_package(Eigen3 3.1.0) 
if (EIGEN3_FOUND) 
	include( ${EIGEN3_USE_FILE} ) 
endif() 


# Add source files in root directory
add_library(${LIB_NAME} SHARED source/CGAL_interface.cpp)

# Linking
set(library_dependencies ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
target_link_libraries(${LIB_NAME} ${library_dependencies})

# Installing
if (NOT ${IS_FORTUNA})
	install (TARGETS ${LIB_NAME} DESTINATION /usr/local/lib)
	install (FILES ${PROJECT_SOURCE_DIR}/include/CGAL_interface.hpp DESTINATION /usr/local/include/CGAL_interface/)
	install (FILES ${PROJECT_SOURCE_DIR}/cmake/CGAL_interfaceConfig.cmake DESTINATION /usr/local/share/CGAL_interface/)
else()
	install (TARGETS ${LIB_NAME} DESTINATION /home/bebe0705/libs/local/lib)
	install (FILES ${PROJECT_SOURCE_DIR}/include/CGAL_interface.hpp DESTINATION /home/bebe0705/libs/local/include/CGAL_interface/)
	install (FILES ${PROJECT_SOURCE_DIR}/cmake/CGAL_interfaceConfig.cmake DESTINATION /home/bebe0705/libs/local/lib/cmake/CGAL_interface/)
endif()

