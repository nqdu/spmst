CXX = g++
OPT= -O3 -march=native -DEIGEN_NO_DEBUG -ffast-math
CXXFLAGS = -std=c++14 -Drestrict=__restrict -Wall 
INCLUDE_DIR = -I./include -I./src/
EIGEN_INC= -I/mnt/c/linuxapp/eigen-3.4.0
PARALLEL = -DEIGEN_DONT_PARALLELIZE -fopenmp

# src and obj
src= $(shell find src/spm2d/*.cpp src/clsqr2/*.cpp src/shared/*.cpp src/spmst/*.cpp )
obj=$(src:%.cpp=%.o)

out=./bin/tomo

all: $(out)

$(out):$(obj)
	$(CXX) -o $(out) $(obj) $(LIB_DIR) $(OPT) $(PARALLEL)

%.o:%.cpp
	$(CXX) -g -c $(CXXFLAGS) $(EIGEN_INC) $(INCLUDE_DIR) $< -o $@ $(OPT) $(PARALLEL)

clean:
	rm -f src/spm2d/*.o src/clsqr2/*.o src/shared/*.o src/spmst/*.o
