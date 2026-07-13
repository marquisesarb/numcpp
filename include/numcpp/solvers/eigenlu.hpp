#pragma once 
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace numcpp::solvers {

    inline Eigen::VectorXd eigenMatrixSolverLU(const Eigen::MatrixXd& A, const Eigen::VectorXd& d) {

        return A.lu().solve(d);
    }



}