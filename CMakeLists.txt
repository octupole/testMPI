# Policy
cmake_minimum_required( VERSION 3.0) 


# Project name
project(ForceField CXX)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -fno-math-errno")  
option(PARALLEL "Enable Parallel compilation " ON)

if (PARALLEL)
   find_package(MPI)
   if (NOT MPI_CXX_FOUND)
 	  message(STATUS "========== OpenTRAJ Parallel Information ==========")
	  message(STATUS )
      message(STATUS "No MPI compiler found. Perhaps you want to set one explicitly?")
      message(STATUS "To override the default compiler, set the environment variable")
      message(STATUS "     export CXX=/full/path/to/mpic++")
      message(STATUS "in bash or ksh or sh.  In csh or tcsh, use:")
      message(STATUS "     setenv CXX /full/path/to/mpic++")
      message(STATUS "before running the cmake command.")
      message(STATUS )
   	  message(STATUS "==================================================")
   endif()
   IF ( MPI_CXX_FOUND )
   	  if(${CMAKE_MINOR_VERSION} GREATER_EQUAL 10)
   	  	if(NOT MPI_MPICXX_FOUND)
   	  		message(FATAL_ERROR "Your MPI implementation does not accept MPI c++ bindings. Connot compile.")
   	  	endif()
   	  elseif()
   	  	message(WARNING "Cannot figure out if your MPI implementation accept MPI c++ bindings. If it does not, it will not compile. Wait and see.")
   	  endif()
      INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
      if(DEFINED $ENV{CXX})
        message(STATUS "Found MPI compiler, but using the CXX compiler set in the CXX environment variable, $ENV{CXX}")
      else()
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
        set(ENV{CXX} ${MPI_CXX_COMPILER})
      endif()
	  add_definitions(-DHAVE_MPI=1)
   ENDIF()
   include(FindOpenMP)
endif()

find_package(Boost 1.35 REQUIRED COMPONENTS
             mpi serialization)

include_directories(${Boost_INCLUDE_DIRS})
include_directories("Parallel")

add_subdirectory(Parallel)

add_executable(tstMPI.x main.cpp)

target_link_libraries(ffield.x InputIO forceField Sys MyTools parallel ${Boost_LIBRARIES})
message( STATUS "LIBRARIES ........... = ${Boost_LIBRARIES}")
message( STATUS "DIRS      ........... = ${Boost_INCLUDE_DIRS}")
