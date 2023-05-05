#!/bin/env bash

work_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

STENCILS_BUILDS=(
    '["XTENSOR_1D","ARRAY_1D","EIGEN_1D","STD_VECTOR_1D","EIGEN_VECTOR"]'
    '["XTENSOR_2D","ARRAY_2D","EIGEN_2D","STD_VECTOR_2D","ARMADILLO","BOOST_UBLAS","BLAZE"]'
    '["VEXCL_CPU","VEXCL_GPU","VIENNACL_WITH_OPENMP","VIENNACL_WITH_OPENCL"]'
)

export FILE_TYPES='["eps", "png", "svg"]'

eval "$(conda shell.bash hook)" && \
conda activate bench-env

results_nbr=1
for stencils_build in ${STENCILS_BUILDS[@]};
do
    export STENCILS_BUILD=$stencils_build
    python3 $PWD/read_results.py "stencil" $((results_nbr++))
done