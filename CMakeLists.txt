#  CMakeLists.txt

cmake_minimum_required(VERSION 3.10.0)

################## compiler options ##################
# compiler
set(CMAKE_CXX_COMPILER          ${CXX})
set(CMAKE_FORTRAN_COMPILER     ${FC})

# compiler flags
set(CMAKE_CXX_FLAGS             ${CXXFLAGS}) 

################## dependencies ##################
# preferred installation prefix of dependencies
set(EIGEN_INC                 ${EIGEN_INC})

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()
setDefault(EIGEN_INC   "/Users/nqdu/software/eigen-3.4.0")

# project
project(SPMTOMO LANGUAGES CXX Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/obj)

################# cmake setup #################
set(CMAKE_CXX_FLAGS  "-Wall -g -DEIGEN_DONT_PARALLELIZE")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") 
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native -ffast-math")
  #set(CMAKE_CXX_FLAGS_RELEASE  "-O3")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -xHost -DEIGEN_STRONG_INLINE=inline")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -extend-source")
endif()

set(CMAKE_VERBOSE_MAKEFILE off)

################# find packages #################
# path of *.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/include/cmake)

#find openmp
find_package(OpenMP)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

# eigen
include_directories(${EIGEN_INC})

################# include paths #################
include_directories(${PROJECT_SOURCE_DIR}/src  ${PROJECT_SOURCE_DIR}/include )

# subs
add_subdirectory(${PROJECT_SOURCE_DIR}/src/shared)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/clsqr2)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/SWD)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/spm2d)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/spmst3D)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/spmst2D)

# exe
add_executable(tomo3D ${PROJECT_SOURCE_DIR}/src/spmst3D/main.cpp)
add_executable(syn ${PROJECT_SOURCE_DIR}/src/spm2d/main.cpp)
add_executable(tomo2D ${PROJECT_SOURCE_DIR}/src/spmst2D/main.cpp)
target_link_libraries(tomo3D spmst3D spm2d swd share lsqr OpenMP::OpenMP_CXX)
target_link_libraries(syn spm2d share lsqr OpenMP::OpenMP_CXX)
target_link_libraries(tomo2D spmst2D spm2d swd share lsqr OpenMP::OpenMP_CXX)