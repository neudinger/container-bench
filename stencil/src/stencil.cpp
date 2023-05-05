#include "stencil.hpp"
// numcpp
void run(std::size_t DIM_SIZE, std::size_t MAX_IT)
{
    const std::string params = std::to_string(DIM_SIZE) + "-" + std::to_string(MAX_IT);
// =====================================================
// VexCL
// GPU / CPU
#if defined(RUN_ARRAY_1D)
    auto marker = ("ARRAY_1D-" + params);
    // std::cout << marker << std::endl;
    ARRAY_1D(marker.c_str());
    std::cout << matrix_tmp[0] << std::endl;
    // for (size_t i = 0; i < DIM_SIZE; ++i)
    // {
    //     for (size_t j = 0; j < DIM_SIZE; ++j)
    //         std::cout << matrixtmp[j + i * DIM_SIZE] << "\t";
    //     std::cout << std::endl;
    // }
#endif // RUN_ARRAY_1D
#if defined(RUN_ARRAY_2D)
    auto marker = ("ARRAY_2D-" + params);
    // std::cout << marker << std::endl;
    ARRAY_2D(marker.c_str());
    // for (auto rowIdx = 0; rowIdx != DIM_SIZE; ++rowIdx)
    // {
    //     for (auto colIdx = 0; colIdx != DIM_SIZE; ++colIdx)
    //         std::cout << matrix_tmp[rowIdx][colIdx] << "\t";
    //     std::cout << std::endl;
    // }
#endif // RUN_ARRAY_2D
#if defined(RUN_BLAZE)
    auto marker = ("BLAZE-" + params);
    // std::cout << marker << std::endl;
    BLAZE(marker.c_str());
#endif // RUN_BLAZE
#if defined(RUN_VEXCL_CPU)
    auto marker = ("RUN_VEXCL_CPU-" + params);
    VexCL(CPU, marker.c_str());
#endif // RUN_VEXCL_CPU
#if defined(RUN_VEXCL_GPU)
    auto marker = ("RUN_VEXCL_GPU-" + params);
    VexCL(GPU, marker.c_str());
#endif // RUN_VEXCL_GPU
// =====================================================
// Viennacl
#if defined(RUN_VIENNACL_WITH_OPENCL)
    auto marker = ("RUN_VIENNACL_WITH_OPENCL-" + params);
    // std::cout << "VIENNACL_WITH_OPENCL " + params << std::endl;
    ViennaCL(marker.c_str());
#endif // VIENNACL_WITH_OPENCL
#if defined(RUN_VIENNACL_WITH_OPENMP)
    auto marker = ("RUN_VIENNACL_WITH_OPENMP-" + params);
    // std::cout << "VIENNACL_WITH_OPENMP " + params << std::endl;
    ViennaCL(marker.c_str());
#endif // VIENNACL_WITH_OPENMP
// =====================================================
// STD::VECTOR
#if defined(RUN_STD_VECTOR_1D)
    auto marker = ("VECTOR_1D-" + params);
    // std::cout << marker << std::endl;
    VECTOR_1D(marker.c_str());
    // for (size_t i = 0; i < DIM_SIZE; ++i)
    // {
    //     for (size_t j = 0; j < DIM_SIZE; ++j)
    //         std::cout << matrix[j + i * DIM_SIZE] << "\t";
    //     std::cout << std::endl;
    // }
#endif // VECTOR1D
#if defined(RUN_STD_VECTOR_2D)
    auto marker = ("VECTOR_2D-" + params);
    // std::cout << marker << std::endl;
    VECTOR_2D(marker.c_str());
#endif // VECTOR2D
       // =====================================================
       // Eigen
#if defined(RUN_EIGEN_1D)
    auto marker = ("Eigen_1D-" + params);
    // std::cout << marker << std::endl;
    EIGEN_1D(marker.c_str());
#endif // EIGEN_1D
#if defined(RUN_EIGEN_2D)
    auto marker = ("Eigen_2D-" + params);
    // std::cout << marker << std::endl;
    EIGEN_2D(marker.c_str());
#endif // EIGEN
#if defined(RUN_EIGEN_VECTOR)
    auto marker = ("EIGEN_VECTOR-" + params);
    // std::cout << marker << std::endl;
    EIGEN_VECTOR(marker.c_str());
    std::cout << matrix_tmp[DIM_SIZE + 1] << std::endl;
#endif // RUN_EIGEN_VECTOR
       // =====================================================
       // BOOST ublas
#if defined(RUN_BOOST_UBLAS)
    auto marker = ("BOOST_ublas-" + params);
    // std::cout << marker << std::endl;
    BOOST_UBLAS(marker.c_str());
    // for (unsigned i = 0; i < system_boost.size1(); ++i)
    //     std::cout << boostblas::row(system_boost, i) << std::endl;
#endif // RUN_BOOST_UBLAS
// =====================================================
// Armadillo
#if defined(RUN_ARMADILLO)
    auto marker = ("Armadillo-" + params);
    // std::cout << marker << std::endl;
    ARMADILLO(marker.c_str());
#endif // RUN_ARMADILLO
// =====================================================
// XTENSOR
#if defined(RUN_XTENSOR_1D)
    auto marker = ("Xtensor_1D-" + params);
    // std::cout << marker << std::endl;
    XTENSOR_1D(marker.c_str());
#endif // RUN_XTENSOR1D
#if defined(RUN_XTENSOR_2D)
    auto marker = ("Xtensor_2D-" + params);
    // std::cout << marker << std::endl;
    XTENSOR_2D(marker.c_str());
#endif // RUN_XTENSOR
    // std::cout << matrix_tmp << std::endl;
}

int main(int argc, char **argv)
{
    ssize_t DIM_SIZE = 10, MAX_IT = 1;
    if (argc >= 3)
    {
        DIM_SIZE = atoll(argv[1]);
        MAX_IT = atoll(argv[2]);

        if (!(DIM_SIZE && MAX_IT))
        {
            std::cerr << "Input parameters errors" << std::endl;
            exit(-1);
        }
    }
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    run(DIM_SIZE, MAX_IT);
    LIKWID_MARKER_CLOSE;
    return 0;
}