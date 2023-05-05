# Bench of c++ popular containers array with likwid

![](https://img.shields.io/badge/Ubuntu-E95420?style=for-the-badge&logo=ubuntu&logoColor=white) ![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)

![C++](https://img.shields.io/badge/c++17-%2300599C.svg?style=for-the-badge&logo=c%2B%2B&logoColor=white) ![](https://img.shields.io/badge/Python_3-FFD43B?style=for-the-badge&logo=python&logoColor=blue) ![](https://img.shields.io/badge/Shell_Script-121011?style=for-the-badge&logo=gnu-bash&logoColor=white)

![](https://img.shields.io/badge/CMake-064F8C?style=for-the-badge&logo=cmake&logoColor=white) ![](https://img.shields.io/badge/conda_env-342B029.svg?&style=for-the-badge&logo=anaconda&logoColor=white)


![](https://img.shields.io/badge/Plotly-239120?style=for-the-badge&logo=plotly&logoColor=white) ![](https://img.shields.io/badge/Pandas-2C2D72?style=for-the-badge&logo=pandas&logoColor=white) ![](https://img.shields.io/badge/Numpy-777BB4?style=for-the-badge&logo=numpy&logoColor=white)

![](https://hpc.fau.de/files/2017/10/ll2-285x200.png) ![](https://upload.wikimedia.org/wikipedia/commons/c/cd/Boost.png) ![](https://mamba.readthedocs.io/en/latest/_static/logo.png)


## Install commands

Required


- Opencl header in `/usr/lib64/`
- C++ Compiler gnu g++ or LLVM clang++ 
- [micromamba](https://mamba.readthedocs.io/en/latest/installation.html#micromamba)

```bash
mamba env create -f environment.yml
conda activate bench-env

export C_INCLUDE_PATH=${C_INCLUDE_PATH}:${CONDA_PREFIX}/include
export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:${CONDA_PREFIX}/include
```

- Download and install the required dependencies 
  ```bash
  ./dl_lib.sh
  ```

- Build every binary
  ```bash
  ./build.sh
  ```

- Run every binary in ./binaries directory with likwid monitoring and write csv result in ./results/csv/
  ```bash
  ./run.sh
  ```

- Read csv to convert in eps, png, svg image performance chart
  ```
  ./read_results.sh
  ```

## Likwid

```sh
# https://rrze-hpc.github.io/likwid/Doxygen/likwid-perfctr.html
likwid-perfctr -O -m -C S0:0 -g MEM_DP ./binaries/stencil_EIGEN_1D 10 1000 > EIGEN_1D.csv

likwid-perfctr -O -m -C S0:0 -g MEM_DP ./binaries/map_STD-3RK\=int64_t-MAP_VALUE\=k1 10 > map_STD\=int64_t-MAP_VALUE\=k1.csv
```

## Stencil

The stencil operation is 2D jacobi stencil

See some result [result](./result.md)

Bench between most popular c++ linear algebra library

- [x] VEXCL
  - [x] CPU
  - [x] GPU
- [x] VIENNACL
  - [x] openCL
  - [x] openMP
- [x] STD VECTOR
  - [x] 1D
  - [x] 2D
- [x] EIGEN
  - [x] 1D
  - [x] 2D
  - [x] VECTOR
- [x] ARMADILLO
- [x] XTENSOR
  - [x] 1D
  - [x] 2D
- [x] BOOST UBLAS
- [x] BLAZE
- [x] ARRAY (double)
  - [x] 1D
  - [x] 2D

<!-- https://img.shields.io/badge/OpenGL-FFFFFF?style=for-the-badge&logo=opengl -->
<!-- - [ ] native GPU (incoming)
  - [ ] Thrust
  - [ ] CUDA
  - [ ] OpenCL
- [ ] Optimised CPU (incoming)
 - [ ] OpenMP -->