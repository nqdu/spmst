CXX = g++
OPT= -O3 -march=native -DEIGEN_NO_DEBUG -ffast-math
CXXFLAGS = -std=c++14 -Drestrict=__restrict -Wall 
INCLUDE_DIR = -I./include -I./src/
EIGEN_INC= -I/mnt/c/linuxapp/eigen-3.4.0
PARALLEL = -DEIGEN_DONT_PARALLELIZE -fopenmp

# src and obj
srcbase= $(shell find src/spm2d/*.cpp|grep -v main.cpp) $(shell find src/clsqr2/*.cpp)  \
	 $(shell find src/shared/*.cpp) $(shell find src/spmst/*.cpp |grep -v main.cpp)
src1 = $(srcbase) src/spmst/main.cpp
src2 = $(srcbase) src/spm2d/main.cpp
obj1 = $(src1:%.cpp=%.o)
obj2 = $(src2:%.cpp=%.o)
prog1 = ./bin/tomo
prog2 = ./bin/travel

all: $(prog1) $(prog2)

$(prog1):$(obj1)
	$(CXX) -o $(prog1)  $(obj1) $(LIB_DIR) $(OPT) $(PARALLEL)

$(prog2):$(obj2)
	$(CXX) -o $(prog2)  $(obj2) $(LIB_DIR) $(OPT) $(PARALLEL)

%.o:%.cpp
	$(CXX) -g -c $(CXXFLAGS) $(EIGEN_INC) $(INCLUDE_DIR) $< -o $@ $(OPT) $(PARALLEL)

clean:
	rm -f src/spm2d/*.o src/clsqr2/*.o src/shared/*.o src/spmst/*.o
