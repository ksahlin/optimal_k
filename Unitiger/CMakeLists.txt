project(Unitiger)

cmake_minimum_required(VERSION 2.6)

################################################################################
# Check for lambda expressions
################################################################################
INCLUDE(CheckCXXSourceCompiles)

IF (CMAKE_COMPILER_IS_GNUXX)
  SET(CMAKE_REQUIRED_FLAGS "-std=c++0x")
ENDIF(CMAKE_COMPILER_IS_GNUXX)

#CHECK_CXX_SOURCE_COMPILES(
#"int main(int argc, char *args[])
#{
#int ret=0;
#return ([&ret]() -> int { return ret; })();
#}
#" CXX_HAS_LAMBDA)

#IF (NOT CXX_HAS_LAMBDA)
#  MESSAGE (FATAL_ERROR "Lambda expressions not available. Use a newer C++ compiler (e.g. GCC Version 4.5 or greater)")
#ENDIF (NOT CXX_HAS_LAMBDA)

################################################################################
# Define cmake modules directory
################################################################################
FOREACH (path "cmake" "../cmake"  "thirdparty/gatb-core/cmake"  "../../thirdparty/gatb-core/gatb-core/cmake")
IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/")
SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
ENDIF()
ENDFOREACH(path)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
#SET (GATB_CORE_EXCLUDE_TOOLS     1)
#SET (GATB_CORE_EXCLUDE_TESTS     1)
#SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)


################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# we give the headers directories from : 
#       - from project source 
#       - from GATB-CORE source
#       - from dsk source 
include_directories (src  ${gatb-core-includes})

# we define the files to be compiled
file (GLOB_RECURSE  ProjectFiles  src/*)

# we define the artifact to be built: the project binary
add_executable        (${CMAKE_PROJECT_NAME}  ${ProjectFiles})

# we define which libraries to be linked with project binary
target_link_libraries (${CMAKE_PROJECT_NAME}  ${gatb-core-libraries})

# we add a compiler flag
SET(GCC_COMPILE_MAP_FLAGS "-std=c++0x -g -fopenmp")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_MAP_FLAGS}")
