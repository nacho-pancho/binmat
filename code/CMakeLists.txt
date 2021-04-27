cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_COMPILER "g++")
project (BINMAT)

set (BINMAT_VERSION_MAJOR 1)
set (BINMAT_VERSION_MINOR 0)

option (PARALLEL "Use parallel instructions." OFF)

if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
  message(STATUS "Target is 64 bits")
  set(TARGETBITS "64")
else("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
  message(STATUS "Target is 32 bits")
  set(TARGETBITS "32")
endif("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")


set (CMAKE_FLAGS   "-mtune=native -march=native -Wformat -Wparentheses -fPIC")
set (CMAKE_CXX_FLAGS " ${CMAKE_FLAGS}")
if(PARALLEL)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(EXTLIB "${EXTLIB} -lgomp")
  add_definitions(-DPARALLEL)
endif()
set(EXTLIB "-lm -lgsl")
set (CMAKE_CXX_FLAGS_DEBUG     "${CMAKE_CXX_FLAGS} -Wall -Wextra -DDEBUG -g -O0 -fno-omit-frame-pointer"  )
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

file(GLOB LIBSRC . lib/bm_*.cpp)

add_library(binmat ${LIBSRC})

include_directories(
  ${PROJECT_SOURCE_DIR}/lib
)

add_custom_command(
  TARGET binmat 
  POST_BUILD 
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/libbinmat.a ${PROJECT_SOURCE_DIR}/lib)

include_directories( ${PROJECT_SOURCE_DIR}/binmat )

add_subdirectory(tools)
