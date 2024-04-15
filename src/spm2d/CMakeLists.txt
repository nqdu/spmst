if(${USE_CUDA})
    file(GLOB SOURCE "*.cpp" "*.cu")
else()
    file(GLOB SOURCE "*.cpp")
endif()

#message("myfile = ${SOURCE}")
#file(REMOVE ${SOURCE} main.cpp )
list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

#message("myfile = ${SOURCE}")

add_library(spm2d ${SOURCE})