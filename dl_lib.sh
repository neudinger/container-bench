#!/bin/env bash

declare DEV_WORKDIR=${PWD}
declare libdir=${PWD}/libs

declare ALL_LIBS=(
    https://github.com/ddemidov/vexcl/archive/refs/tags/1.4.2.tar.gz
    https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    http://sourceforge.net/projects/arma/files/armadillo-10.7.0.tar.xz
    http://sourceforge.net/projects/viennacl/files/1.7.x/ViennaCL-1.7.1.tar.gz/download
    https://github.com/KhronosGroup/OpenCL-Headers/archive/refs/tags/v2021.06.30.tar.gz
    https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.bz2
    https://bitbucket.org/blaze-lib/blaze/downloads/blaze-3.8.tar.gz
)

mkdir -p ${libdir};

for lib in ${ALL_LIBS[@]};
do
    curl --progress-bar -L ${lib} -o lib.tar.gz && \
    tar -xvf lib.tar.gz -C ${libdir} && \
    rm lib.tar.gz
done

# # ====BUILD LIB====
export PREFIX=${CONDA_PREFIX}
export CMAKE_INSTALL_PREFIX=${CONDA_PREFIX}

# Install boost
cd ${libdir}/boost_1_77_0/
./bootstrap.sh
./b2


# Install LIKWID
cd ${DEV_WORKDIR};

# Install vexcl
cmake -S ${libdir}/vexcl-1.4.2 -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
-DCMAKE_BUILD_TYPE=Release \
-DBoost_NO_BOOST_CMAKE:BOOL=TRUE \
-DBoost_NO_SYSTEM_PATHS:BOOL=TRUE \
-DBOOST_ROOT=${libdir}/boost_1_77_0 \
-DOpenCL_FOUND=True -DOpenCL_INCLUDE_DIR=/usr/lib64/ -DUSE_LIBCPP:BOOL=TRUE \
-B ${libdir}/vexcl-1.4.2/build ${libdir}/vexcl-1.4.2

cd ${libdir}/vexcl-1.4.2/build && make install;


# Install LIKWID
cd ${DEV_WORKDIR};
# sed in file because with sudo make install PREFIX env disappear
wget https://github.com/RRZE-HPC/likwid/archive/refs/tags/v5.2.0.tar.gz && \
tar -xvf v5.2.0.tar.gz -C ${libdir} && rm v5.2.0.tar.gz && \
sed -i "s|\/usr\/local|$CONDA_PREFIX|" ${libdir}/likwid-5.2.0/config.mk && \
cd ${libdir}/likwid-5.2.0 && make && sudo make install

cd ${DEV_WORKDIR};
