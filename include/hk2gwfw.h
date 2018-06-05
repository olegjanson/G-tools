#include <iostream>
#include <fstream>
#include <vector>
#if __has_include("mklaa.h")
#define EIGEN_USE_MKL_ALL
#endif
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
