#  CMakeLists.txt

cmake_minimum_required(VERSION 3.20.0)

################## compiler options ##################
# compiler
set(CMAKE_CXX_COMPILER          ${CXX})
set(CMAKE_FORTRAN_COMPILER     ${FC})

# compiler flags
set(CMAKE_CXX_FLAGS             ${CXXFLAGS}) 

################## dependencies ##################
# preferred installation prefix of dependencies
set(EIGEN_INC                 ${EIGEN_INC})
set(SPMST_NPT_AUX               ${SPMST_NPT_AUX})

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()
setDefault(EIGEN_INC   "~/software/eigen-5.0.0")
setDefault(SPMST_NPT_AUX 9)

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
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  #set(CMAKE_CXX_FLAGS_RELEASE  "-O3")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -xHost -DEIGEN_STRONG_INLINE=inline")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -extend-source")
else()
  # this is macos ,use same as gnu
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
endif()

# add define macro for NPT_AUX
add_definitions(-DSPMST_NPT_AUX=${SPMST_NPT_AUX})

set(CMAKE_VERBOSE_MAKEFILE off)

################# find packages #################
# path of *.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/include/cmake)

#find openmp
find_package(OpenMP)
if (OPENMP_FOUND)
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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
add_executable(tomo3d ${PROJECT_SOURCE_DIR}/src/spmst3D/main.cpp)
add_executable(syn2d ${PROJECT_SOURCE_DIR}/src/spm2d/main.cpp)
add_executable(tomo2d ${PROJECT_SOURCE_DIR}/src/spmst2D/main.cpp)
target_link_libraries(tomo3d spmst3D spm2d swd share lsqr OpenMP::OpenMP_CXX)
target_link_libraries(syn2d spm2d share lsqr OpenMP::OpenMP_CXX)
target_link_libraries(tomo2d spmst2D spm2d swd share lsqr OpenMP::OpenMP_CXX)