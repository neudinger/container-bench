#if !defined(STENCIL)
#define STENCIL

#include <iostream>

#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

// https://man7.org/linux/man-pages/man3/alloca.3.html
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#ifdef _WIN32
// windows code goes here
#include <windows.h>
#define malloca _alloca
#else
// posix code goes here
#include <alloca.h>
#define malloca alloca
#endif

#if defined(RUN_VECTOR)
#include <vector>
#endif

#define CL_TARGET_OPENCL_VERSION 300

#if defined(RUN_BLAZE)
#include <blaze/Blaze.h>
#include <blaze/Forward.h>
#endif

#if defined(RUN_VEXCL)
#include <boost/compute.hpp>
#include <vector>
#include <vexcl/vexcl.hpp>
#endif // RUN_VEXCL

#if defined(RUN_VIENNACL)
#include <viennacl/matrix.hpp>
#include <viennacl/matrix_proxy.hpp>
#endif // RUN_VIENNACL

#if defined(RUN_ARMADILLO)
#include <armadillo>
#endif // RUN_ARMADILLO

#if defined(RUN_EIGEN)
#include <Eigen/Core>
#include <Eigen/Dense>
#endif // RUN_EIGEN

#if defined(RUN_XTENSOR)
#include <xtensor/xio.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xview.hpp>
#endif // RUN_XTENSOR

#if defined(RUN_BOOST_UBLAS)
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace boostblas = boost::numeric::ublas;
#endif // RUN_BOOST_UBLAS

#define xstr(s) str(s)
#define str(s) #s

#define ARRAY_1D(marker)                                                                                               \
    double *matrix = (double *)malloca((DIM_SIZE * DIM_SIZE) * sizeof(double));                                        \
    double *matrix_tmp = (double *)malloca((DIM_SIZE * DIM_SIZE) * sizeof(double));                                    \
    for (size_t i = 0; i < DIM_SIZE; i++)                                                                              \
    {                                                                                                                  \
        matrix[i] = 100;                                                                                               \
        matrix_tmp[i] = 100;                                                                                           \
    }                                                                                                                  \
    for (size_t i = DIM_SIZE; i < DIM_SIZE * DIM_SIZE; i++)                                                            \
    {                                                                                                                  \
        matrix[i] = 0;                                                                                                 \
        matrix_tmp[i] = 0;                                                                                             \
    }                                                                                                                  \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        for (size_t i = 1; i < DIM_SIZE - 1; ++i)                                                                      \
            for (size_t j = 1; j < DIM_SIZE - 1; ++j)                                                                  \
                matrix_tmp[j + i * DIM_SIZE] =                                                                         \
                    0.25 * (matrix[j + i * DIM_SIZE + 1] + matrix[j + i * DIM_SIZE - 1] +                              \
                            matrix[j + i * DIM_SIZE - DIM_SIZE] + matrix[j + i * DIM_SIZE + DIM_SIZE]);                \
    LIKWID_MARKER_STOP(marker);

#define ARRAY_2D(marker)                                                                                               \
    double *ptr, *ptr_tmp, **matrix, **matrix_tmp;                                                                     \
    ssize_t len = sizeof(double *) * DIM_SIZE + sizeof(double) * DIM_SIZE * DIM_SIZE;                                  \
    matrix = (double **)malloca(len);                                                                                  \
    matrix_tmp = (double **)malloca(len);                                                                              \
    ptr = (double *)(matrix + DIM_SIZE);                                                                               \
    ptr_tmp = (double *)(matrix_tmp + DIM_SIZE);                                                                       \
    for (size_t i = 0; i < DIM_SIZE; i++)                                                                              \
    {                                                                                                                  \
        matrix[i] = (ptr + DIM_SIZE * i);                                                                              \
        matrix_tmp[i] = (ptr_tmp + DIM_SIZE * i);                                                                      \
    }                                                                                                                  \
    for (size_t i = 0; i < DIM_SIZE; i++)                                                                              \
    {                                                                                                                  \
        matrix[0][i] = 100;                                                                                            \
        matrix_tmp[0][i] = 100;                                                                                        \
    }                                                                                                                  \
    for (size_t rowIdx = 1; rowIdx != DIM_SIZE; ++rowIdx)                                                              \
        for (size_t colIdx = 0; colIdx != DIM_SIZE; ++colIdx)                                                          \
        {                                                                                                              \
            matrix[rowIdx][colIdx] = 0;                                                                                \
            matrix_tmp[rowIdx][colIdx] = 0;                                                                            \
        }                                                                                                              \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        for (size_t rowIdx = 1; rowIdx != DIM_SIZE - 1; ++rowIdx)                                                      \
            for (size_t colIdx = 1; colIdx != DIM_SIZE - 1; ++colIdx)                                                  \
                matrix_tmp[rowIdx][colIdx] = 0.25 * (matrix[rowIdx - 1][colIdx] + matrix[rowIdx + 1][colIdx] +         \
                                                     matrix[rowIdx][colIdx - 1] + matrix[rowIdx][colIdx + 1]);         \
    LIKWID_MARKER_STOP(marker);

#define VexCL(TYPE, marker)                                                                                            \
    std::vector<double> matrix_stdVector(DIM_SIZE *DIM_SIZE, 0);                                                       \
    for (size_t i = 0; i < DIM_SIZE; i++)                                                                              \
        matrix_stdVector[i] = 100;                                                                                     \
    std::cout << "VexCL " str(TYPE) " : " << std::endl;                                                                \
    vex::Context ctx(vex::Filter::TYPE &&vex::Filter::DoublePrecision);                                                \
    if (!ctx)                                                                                                          \
        throw std::runtime_error("No devices available.");                                                             \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    vex::backend::kernel kernel_apply(                                                                                 \
        ctx.queue(0),                                                                                                  \
        VEX_STRINGIZE_SOURCE(kernel void kernel_apply(ulong _dim, ulong _iteration, global double *_system,            \
                                                      global double *_system_tmp) {                                    \
            size_t iloc = get_local_id(0);                                                                             \
            size_t nloc = get_local_size(0);                                                                           \
            for (ulong max_iteration = _iteration; max_iteration > 0; --max_iteration)                                 \
                for (size_t i = iloc + 1; i < _dim - 1; i += nloc)                                                     \
                    for (size_t j = 1; j < _dim - 1; ++j)                                                              \
                        _system_tmp[j + i * _dim] =                                                                    \
                            0.25 * (_system[j + i * _dim + 1] + _system[j + i * _dim - 1] +                            \
                                    _system[j + i * _dim - _dim] + _system[j + i * _dim + _dim]);                      \
        }),                                                                                                            \
        "kernel_apply");                                                                                               \
    vex::vector<double> matrix(ctx, matrix_stdVector);                                                                 \
    vex::vector<double> matrix_tmp(ctx, matrix_stdVector);                                                             \
    kernel_apply(ctx.queue(0), static_cast<cl_ulong>(DIM_SIZE), static_cast<cl_ulong>(MAX_IT), matrix(0),              \
                 matrix_tmp(0));                                                                                       \
    LIKWID_MARKER_STOP(marker);

#define ViennaCL(marker)                                                                                               \
    viennacl::matrix<double> matrix(DIM_SIZE, DIM_SIZE);                                                               \
    viennacl::matrix<double> matrix_tmp(DIM_SIZE, DIM_SIZE);                                                           \
    matrix.clear();                                                                                                    \
    for (std::size_t j = 0; j < DIM_SIZE; ++j)                                                                         \
        matrix(0, j) = 100;                                                                                            \
    const std::size_t row_size = matrix.size1();                                                                       \
    const std::size_t cols_size = matrix.size2();                                                                      \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        viennacl::matrix_range<viennacl::matrix<double>>(matrix_tmp, viennacl::range(1, row_size - 1),                 \
                                                         viennacl::range(1, cols_size - 1)) =                          \
            0.25 * (viennacl::matrix_range<viennacl::matrix<double>>(matrix, viennacl::range(0, row_size - 2),         \
                                                                     viennacl::range(1, cols_size - 1)) +              \
                    viennacl::matrix_range<viennacl::matrix<double>>(matrix, viennacl::range(2, row_size),             \
                                                                     viennacl::range(1, cols_size - 1)) +              \
                    viennacl::matrix_range<viennacl::matrix<double>>(matrix, viennacl::range(1, row_size - 1),         \
                                                                     viennacl::range(0, cols_size - 2)) +              \
                    viennacl::matrix_range<viennacl::matrix<double>>(matrix, viennacl::range(1, row_size - 1),         \
                                                                     viennacl::range(2, cols_size)));                  \
    LIKWID_MARKER_STOP(marker);

#define VECTOR_2D(marker)                                                                                              \
    using MatrixVectorDb = std::vector<std::vector<double>>;                                                           \
    MatrixVectorDb matrix(DIM_SIZE, std::vector<double>(DIM_SIZE, 0));                                                 \
    MatrixVectorDb matrix_tmp(DIM_SIZE, std::vector<double>(DIM_SIZE, 0));                                             \
    for (auto &col : matrix[0])                                                                                        \
        col = 100;                                                                                                     \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        for (auto rowIdx = 1; rowIdx != DIM_SIZE - 1; ++rowIdx)                                                        \
            for (auto colIdx = 1; colIdx != DIM_SIZE - 1; ++colIdx)                                                    \
                matrix_tmp[rowIdx][colIdx] = 0.25 * (matrix[rowIdx - 1][colIdx] + matrix[rowIdx + 1][colIdx] +         \
                                                     matrix[rowIdx][colIdx - 1] + matrix[rowIdx][colIdx + 1]);         \
    LIKWID_MARKER_STOP(marker);

#define VECTOR_1D(marker)                                                                                              \
    std::vector<double> matrix(DIM_SIZE *DIM_SIZE, 0);                                                                 \
    for (size_t i = 0; i < DIM_SIZE; i++)                                                                              \
        matrix[i] = 100;                                                                                               \
    std::vector<double> matrix_tmp(DIM_SIZE *DIM_SIZE, 0);                                                             \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        for (size_t i = 1; i < DIM_SIZE - 1; ++i)                                                                      \
            for (size_t j = 1; j < DIM_SIZE - 1; ++j)                                                                  \
                matrix_tmp[j + i * DIM_SIZE] =                                                                         \
                    0.25 * (matrix[j + i * DIM_SIZE + 1] + matrix[j + i * DIM_SIZE - 1] +                              \
                            matrix[j + i * DIM_SIZE - DIM_SIZE] + matrix[j + i * DIM_SIZE + DIM_SIZE]);                \
    LIKWID_MARKER_STOP(marker);

#define ARMADILLO(marker)                                                                                              \
    arma::mat matrix(DIM_SIZE, DIM_SIZE, arma::fill::value<double>(0.0));                                              \
    arma::mat matrix_tmp(DIM_SIZE, DIM_SIZE, arma::fill::value<double>(0.0));                                          \
    arma::rowvec fill(DIM_SIZE, arma::fill::value<double>(100.0));                                                     \
    matrix.row(0) = fill;                                                                                              \
    const std::size_t row_size = matrix.n_rows - 1;                                                                    \
    const std::size_t cols_size = matrix.n_cols - 1;                                                                   \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        matrix_tmp.rows(1, row_size - 1).cols(1, cols_size - 1) =                                                      \
            0.25 *                                                                                                     \
            (matrix.rows(0, row_size - 2).cols(1, cols_size - 1) + matrix.rows(2, row_size).cols(1, cols_size - 1) +   \
             matrix.rows(1, row_size - 1).cols(0, cols_size - 2) + matrix.rows(1, row_size - 1).cols(2, cols_size));   \
    LIKWID_MARKER_STOP(marker);

#define BLAZE(marker)                                                                                                  \
    blaze::DynamicMatrix<double> matrix(DIM_SIZE, DIM_SIZE, 0);                                                        \
    blaze::DynamicMatrix<double> matrix_tmp(DIM_SIZE, DIM_SIZE, 0);                                                    \
    blaze::row(matrix, 0) = 100.;                                                                                      \
    const u_long row_size = DIM_SIZE - 1;                                                                              \
    const u_long cols_size = DIM_SIZE - 1;                                                                             \
    matrix_tmp = matrix;                                                                                               \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        blaze::submatrix(matrix_tmp, 1, 1, row_size - 1, cols_size - 1) =                                              \
            blaze::noalias(0.25 * (blaze::submatrix(matrix, 0, 1, row_size - 1, cols_size - 1) +                       \
                                   blaze::submatrix(matrix, 2, 1, row_size - 1, cols_size - 1) +                       \
                                   blaze::submatrix(matrix, 1, 0, row_size - 1, cols_size - 1) +                       \
                                   blaze::submatrix(matrix, 1, 2, row_size - 1, cols_size - 1)));                      \
    LIKWID_MARKER_STOP(marker);

#define EIGEN_2D(marker)                                                                                               \
    Eigen::MatrixXd matrix(DIM_SIZE, DIM_SIZE);                                                                        \
    Eigen::MatrixXd matrix_tmp(DIM_SIZE, DIM_SIZE);                                                                    \
    const std::size_t row_size = matrix.rows() - 1;                                                                    \
    const std::size_t cols_size = matrix.cols() - 1;                                                                   \
    matrix(0, Eigen::seq(0, DIM_SIZE)) = Eigen::RowVectorXd(DIM_SIZE).setConstant(100);                                \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        matrix_tmp(Eigen::seq(1, row_size - 1), Eigen::seq(1, cols_size - 1)).noalias() =                              \
            0.25 * (matrix(Eigen::seq(0, row_size - 2), Eigen::seq(1, cols_size - 1)) +                                \
                    matrix(Eigen::seq(2, row_size), Eigen::seq(1, cols_size - 1)) +                                    \
                    matrix(Eigen::seq(1, row_size - 1), Eigen::seq(0, cols_size - 2)) +                                \
                    matrix(Eigen::seq(1, row_size - 1), Eigen::seq(2, cols_size)));                                    \
    LIKWID_MARKER_STOP(marker);

#define EIGEN_1D(marker)                                                                                               \
    Eigen::ArrayXd matrix(DIM_SIZE *DIM_SIZE);                                                                         \
    Eigen::ArrayXd matrix_tmp(DIM_SIZE *DIM_SIZE);                                                                     \
    matrix(Eigen::seq(0, DIM_SIZE)) = 100.0;                                                                           \
    const std::size_t row_size = DIM_SIZE;                                                                             \
    const std::size_t cols_size = DIM_SIZE;                                                                            \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        for (size_t j = 1; j < cols_size - 1; ++j)                                                                     \
            matrix_tmp(Eigen::seq(1 + j * row_size, row_size - 1 + j * row_size)) =                                    \
                0.25 * (matrix(Eigen::seq(j * row_size, row_size - 2 + j * row_size)) +                                \
                        matrix(Eigen::seq(2 + j * row_size, row_size + j * row_size)) +                                \
                        matrix(Eigen::seq(1 + (j - 1) * row_size, row_size - 1 + (j - 1) * row_size)) +                \
                        matrix(Eigen::seq(1 + (j + 1) * row_size, row_size - 1 + (j + 1) * row_size)));                \
    LIKWID_MARKER_STOP(marker);

#define EIGEN_VECTOR(marker)                                                                                           \
    Eigen::RowVectorXd matrix(DIM_SIZE *DIM_SIZE);                                                                     \
    Eigen::RowVectorXd matrix_tmp(DIM_SIZE *DIM_SIZE);                                                                 \
    for (size_t i = 0; i < DIM_SIZE; i++)                                                                              \
        matrix[i] = 100;                                                                                               \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        for (size_t i = 1; i < DIM_SIZE - 1; ++i)                                                                      \
            for (size_t j = 1; j < DIM_SIZE - 1; ++j)                                                                  \
                matrix_tmp[j + i * DIM_SIZE] =                                                                         \
                    0.25 * (matrix[j + i * DIM_SIZE + 1] + matrix[j + i * DIM_SIZE - 1] +                              \
                            matrix[j + i * DIM_SIZE - DIM_SIZE] + matrix[j + i * DIM_SIZE + DIM_SIZE]);                \
    LIKWID_MARKER_STOP(marker);

#define XTENSOR_2D(marker)                                                                                             \
    xt::xarray<double> matrix({DIM_SIZE, DIM_SIZE}, 0);                                                                \
    xt::xarray<double> matrix_tmp({DIM_SIZE, DIM_SIZE}, 0);                                                            \
    xt::view(matrix, 0) = 100.;                                                                                        \
    const std::size_t row_size = matrix.shape(0);                                                                      \
    const std::size_t cols_size = matrix.shape(1);                                                                     \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        xt::noalias(xt::view(matrix_tmp, xt::range(1, row_size - 1), xt::range(1, cols_size - 1))) =                   \
            0.25 * (xt::view(matrix, xt::range(0, row_size - 2), xt::range(1, cols_size - 1)) +                        \
                    xt::view(matrix, xt::range(2, row_size), xt::range(1, cols_size - 1)) +                            \
                    xt::view(matrix, xt::range(1, row_size - 1), xt::range(0, cols_size - 2)) +                        \
                    xt::view(matrix, xt::range(1, row_size - 1), xt::range(2, cols_size)));                            \
    LIKWID_MARKER_STOP(marker);

#define XTENSOR_1D(marker)                                                                                             \
    xt::xtensor<double, 1> matrix({DIM_SIZE * DIM_SIZE}, 0);                                                           \
    xt::xtensor<double, 1> matrix_tmp({DIM_SIZE * DIM_SIZE}, 0);                                                       \
    xt::view(matrix, 0) = 100.;                                                                                        \
    const std::size_t row_size = DIM_SIZE;                                                                             \
    const std::size_t cols_size = DIM_SIZE;                                                                            \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        for (size_t j = 1; j < cols_size - 1; ++j)                                                                     \
            xt::noalias(xt::view(matrix_tmp, xt::range(1 + j * row_size, row_size - 1 + j * row_size))) =              \
                0.25 * (xt::view(matrix, xt::range(j * row_size, row_size - 2 + j * row_size)) +                       \
                        xt::view(matrix, xt::range(2 + j * row_size, row_size + j * row_size)) +                       \
                        xt::view(matrix, xt::range(1 + (j - 1) * row_size, row_size - 1 + (j - 1) * row_size)) +       \
                        xt::view(matrix, xt::range(1 + (j + 1) * row_size, row_size - 1 + (j + 1) * row_size)));       \
    LIKWID_MARKER_STOP(marker);

#define BOOST_UBLAS(marker)                                                                                            \
    boostblas::matrix<double> matrix(DIM_SIZE, DIM_SIZE);                                                              \
    boostblas::matrix<double> matrix_tmp(DIM_SIZE, DIM_SIZE);                                                          \
    matrix.clear();                                                                                                    \
    for (unsigned j = 0; j < matrix.size2(); ++j)                                                                      \
        matrix(0, j) = 100;                                                                                            \
    const std::size_t row_size = matrix.size1();                                                                       \
    const std::size_t cols_size = matrix.size2();                                                                      \
    LIKWID_MARKER_REGISTER(marker);                                                                                    \
    LIKWID_MARKER_START(marker);                                                                                       \
    for (size_t max_iteration = MAX_IT; max_iteration > 0; --max_iteration)                                            \
        boostblas::noalias(boostblas::subrange(matrix_tmp, 1, row_size - 1, 1, cols_size - 1)) =                       \
            0.25 * (boostblas::subrange(matrix, 0, row_size - 2, 1, cols_size - 1) +                                   \
                    boostblas::subrange(matrix, 2, row_size, 1, cols_size - 1) +                                       \
                    boostblas::subrange(matrix, 1, row_size - 1, 0, cols_size - 2) +                                   \
                    boostblas::subrange(matrix, 1, row_size - 1, 2, cols_size));                                       \
    LIKWID_MARKER_STOP(marker);

#endif // STENCIL