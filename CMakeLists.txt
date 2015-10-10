set(CMAKE_CXX_COMPILER mpic++ ) 
set(CMAKE_LINKER mpic++ )
set(CMAKE_C_COMPILER mpicc )
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpedantic -Wno-reorder -Wno-unused-parameter -fstrict-aliasing ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -Wall -Wextra -Wpedantic -Wno-reorder -Wno-unused-parameter -fstrict-aliasing ")

cmake_minimum_required (VERSION 2.6)
project (nanolens)

set(nanolens_VERSION_MAJOR 2)
set(nanolens_VERSION_MINOR 99)
set(nanolens_VERSION_PATCH 5)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/nanolens.hpp.in"
  "${PROJECT_BINARY_DIR}/nanolens.hpp"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

add_executable(nanolens nanolens.cpp nanolens.xml)
add_executable(nanopostprocess nanopostprocess.cpp)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: Debug Release"
       FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(MPI)

if(MPI_FOUND)
  set(LIBS ${LIBS} ${MPI_C_LIBRARIES})
endif(MPI_FOUND)

find_package(Boost
  1.36      # Minimum or EXACT version e.g. 1.36.0
  REQUIRED             # Fail with error if Boost is not found
  mpi serialization # Boost libraries by their canonical name
  )

set(LIBS ${LIBS} ${BOOST_LIBRARIES})


include_directories(${BOOST_INCLUDE_DIRS})
target_link_libraries (nanolens boost_serialization boost_mpi cfitsio)
target_link_libraries (nanopostprocess boost_serialization boost_mpi cfitsio)

#Install nanolens configuration file
add_custom_command(TARGET nanolens POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                       ${CMAKE_SOURCE_DIR}/nanolens.xml $<TARGET_FILE_DIR:nanolens>)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)



