SPMST: Surface Wave Tomography with Topography
==============

# Table of Contents
1. [Introduction](#Introduction)
2. [Installation](#Installation)
3. [File Format](#File Format)
4. [Fourth Example](#fourth-examplehttpwwwfourthexamplecom)


# 1. Introduction

`SPMST` is a package to conduct 2-D surface wave tomography with topography by using shortest path method. In this method, surface wave is propagating along the (curved) surface.

This manual provides guidance on using the `SPMST` package to perform 2D/3D surface wave ray-based tomography of shear wave velocities. These package includes 3 programs:

- Synthesize phase velocity travel time for a given 2-D phase velocity map and topography.
- 2-D phase velocity tomography (from dispersion to phase maps).
- 3-D direct inversion of both phase/group dispersion curves. 

It should be noted that this package has been tested by a limited number of users under similar computational conditions, so it may contain some bugs. I would appreciate any bug reports submitted on the GitHub page.

For further details, please refer to [here](https://academic.oup.com/gji/article/237/2/1235/7633451).
 

# 2. Installation
Requirements:
name | version|  
|:---:|:----:|  
|[GCC](https://gcc.gnu.org/)| >=7.5|
|[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) | >=3.4.0 | 
|[CMAKE](https://cmake.org/) | >=3.12.0|
|[paraview]()||

After you have downloaded the code from Github, you can compile it by using:
```bash
mkdir build;
cd build;
cmake .. -DCXX=g++ -DFC=gfortran -DEIGEN_INC=/path/to/EIGEN 
make -j8
```
Then all the programs will be put under the directory `bin`.

# 3. File Format
All 3 programs are highly dependent on several input files: topography file, source/receivers file, and the 2D/3D velocity file. Before we go through details in each program, we talk a little bit about the file format:

## 3.1 Topography file
One template of the topography file is like:
```
51 51
-5.000000 5.000000
-5.000000 5.000000
3.72008e-44 
1.87493e-42
8.05255e-41
2.94709e-39
9.19108e-38
2.4426e-36
5.53161e-35
1.06749e-33
```
- The first two numbers are number of points in `x/y` dimension, respectively
- The 2nd/3rd line are minimum/maximum coordinates for `x/y` direction (in km or degree). 
- The following lines are the topography at each point (in meters), it should be written as:
```python
for iy in range(ny):
    for ix in range(nx):
      f.write("%f\n" %(topo[iy,ix]))
```
If you want to use spherical coordinates in your project, you should make sure all coordiantes are in degree, and remember the notation duality: `x -> longitude, y -> latitude`

## 3.2 2D velcity file
here is the template:
```code
40 40
89.200000 90.300000
38.900000 39.300000
1
3.000000
3.000000
3.000000
3.000000
3.000000
3.000000
```
- The first two numbers are number of points in `x/y` dimension, respectively
- The 2nd/3rd line are minimum/maximum coordinates for `x/y` direction (in km or degree). 
- The fourth line is a flag that enables the spherical coordinate system; if set to 1, spherical coordinates will be used.
- The following lines are the velocities at each point (in km/s). The order of for loop is same as the topography file.

## 3.3 3D velocity file
here is the template:
```
33 36 27
99.700000 110.200000
25.700000 35.300000
0 0
0.000000 0.400000 0.800000 1.200000 1.600000 2.000000 2.400000 2.800000 3.200000 3.600000 4.000000 4.400000 4.800000 5.200000 5.600000 6.000000 6.400000 6.800000 7.200000 7.600000 8.000000 8.400000 8.800000 9.200000 9.600000 15.000000 35.000000 
2.800000
2.800000
2.800000
2.800000
2.800000
2.800000
2.800000
```
- The first two numbers are number of points in `x/y/z` dimension, respectively
- The 2nd/3rd line are minimum/maximum coordinates for `x/y` direction. 
- The two numbers in fourth line:
    * `Spherical coordinates flag`, same as 2D case
    * `Shift interface flag`. If it is set to 0, the topography will be processed as the thickening/thinning of the first layer. Otherwise the location of all interface will be shifted up/downward but the relative location is not changed. I `strongly recommend` you to set this flag to 0 unless you really know what you are doing.  
- The 5-th line are the depth for each point (in km).
- The following lines are the velocities at each point (in km/s), it should be written as:
```python
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            f.write("%f\n" % (veloc[iz,iy,ix]))
```

## 3.4 2D Phase velocity file
here is the one block (source-receiver pair) in this file:
```
# 89.488540 39.099072 4
89.987624 39.057909 2.45
89.948282 38.922526 2.54
89.306624 39.269530 2.60
89.425898 39.152494 2.676753
```
This file contains the source/receiver location and phase velocity observations at a single period:
- The first line (must start with `#`) contains the x and y coordinates (in km/deg) for source, and the following number is how may receivers related to this source.
- The following lines are receiver coordiantes (in km/deg) and the observed phase velocity (in km/s). 

You should list all the blocks you required in this file. Another point you should note that this program `cannot` handle group velocity observations.

## 3.5 3D Dispersion File
This file contains the period used for each type of surface wave and the source/receiver location and phase velocity observations at a each period, each wave type. Here is the template:
```
19                                # ntRc (followed by periods)
0.2 0.6 1. 1.4 1.8 2.2 2.6 3. 3.4 3.8 4.2 4.6 5.  5.4 5.8 6.2 6.6 7. 7.4
19                                # ntRg (followed by periods)
0.2 0.6 1.  1.4 1.8 2.2 2.6 3.  3.4 3.8 4.2 4.6 5.  5.4 5.8 6.2 6.6 7  7.4
0                                 # ntLc (followed by periods)
0                                 # ntLg
# 108.646000 31.962600 6 0 2 0
105.764100 33.732800 3.082300
106.020600 34.342500 3.114000
104.991700 33.357400 3.048400
106.139500 33.356800 3.080200
107.817000 34.128300 3.125100
106.800200 33.229300 3.057500
```
Obviously this file contains a header and data blocks. The header parameters are:
- The first line contains the size of the period array for Rayleigh wave phase velocity. The if this number is not 0, the following line include the period vector (in s).
- The following several lines contain the same for Rayleigh wave group velocity, Love wave's phase velocity and Love wave's group velocity, respectively.

The rest of this file is the description for each source receiver pair at a given period of a single wave type. The source is started with a `#`:
```code
# x y nsta pid wave-type disp-type
```
- `x y` coordinates, in km/deg
- `nsta` number of stations for this pair
- `pid` is the index of the period in the period array, start from 0.
- `wave-type` Rayleigh wave (=2) or Love wave (=1)
- `disp-type` phase velocity (=0) or group velocity (=1)
For the template, the source line means: this block has `6` receivers, related to `Rayleigh phase velocity` at `0.2`s.

Then the following `nsta` lines contains the receiver coordinates (same units with source) and the observed dispersion (in km/s).

## 3.6 Inverse problem file
```
# SPMST2D/3D inversion Parameters

# # of iterations
NITERS =  6
ITER_CURRENT = 0 # current iteration

#  minimum velocity, maximum velocity in km/s
MIN_VELOC =  2.0
MAX_VELOC = 5.0 

# synthetic test
SYN_TEST = 1                # synthetic flag(0:real data,1:synthetic)
NOISE_LEVEL = 0.1          #  noise level, std value of gaussian noise

# for LSQR-based inversion 
SMOOTH = 2.0  # 2-nd Tikhonov regularization parameter
DAMP = 0.01    # damping parameter for LSQR
NTHREADS =  2  # # of threads used in LSQR solver 
```
This is a self-explanatory file, you can add any comments in it (start with `#`). To modify it, you should leave at least `one space` in the left/right side of `=`.

# Appendix
This package is trying to minimize the frequency-dependent phase velocity related travel time:
```math
t(\omega) = \min_{P\in all ~ path}\int_{P} \frac{\mathrm{d}s}{c(\omega)}
```
For group velocity, it should be:
```math
t'(\omega) = \int_{P} \frac{\mathrm{d}s}{U(\omega)}
```
where the path $P$ is the `same path` as obtained in corresponding phase velocity case.