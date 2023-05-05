#!/bin/env bash

build_dir=${work_dir}/builds/stencil
here=${work_dir}/stencil

ALL_BUILD=(
    ARMADILLO
    ARRAY_1D
    ARRAY_2D
    BLAZE
    BOOST_UBLAS
    EIGEN_1D
    EIGEN_2D
    STD_VECTOR_1D
    STD_VECTOR_2D
    EIGEN_VECTOR
    VEXCL_CPU
    VEXCL_GPU
    VIENNACL_WITH_OPENCL
    VIENNACL_WITH_OPENMP
    XTENSOR_1D
    XTENSOR_2D
)

FLAGS="-DCMAKE_BUILD_TYPE=Release -DOpenCL_FOUND=True -DOpenCL_INCLUDE_DIR=/usr/lib64/"

for item in ${ALL_BUILD[@]};
do
    #Release OR MinSizeRel OR RelWithDebInfo
    cmake -S ${here} \
    ${FLAGS} \
    -DCMAKE_BUILD_TYPE=MinSizeRel \
    -D${item}:BOOL=ON  \
    -B ${build_dir} && cmake --build ${build_dir} && \
    mv ${build_dir}/stencil ${work_dir}/binaries/stencil_${item};
done
