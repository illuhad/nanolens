set(CMAKE_CXX_COMPILER mpic++ ) 
set(CMAKE_LINKER mpic++ )
set(CMAKE_C_COMPILER mpicc )
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpedantic -Wno-reorder -Wno-unused-parameter -fstrict-aliasing ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -Wall -Wextra -Wpedantic -Wno-reorder -Wno-unused-parameter -fstrict-aliasing ")

cmake_minimum_required (VERSION 2.6)
project (nanolens)
add_executable(nanolens nanolens.cpp)


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



