#pragma once 
#include <Eigen/Dense>
#include <limits>

namespace numcpp::stats {

    struct CorrelationMatrix {

        Eigen::MatrixXd matrix;
        
        CorrelationMatrix(double coeff) {

            matrix = Eigen::MatrixXd::Zero(2,2);
            coeff = std::abs(coeff) <= 1.0 ? coeff : std::numeric_limits<double>::quiet_NaN();
            matrix << 1.0, coeff, 
                    coeff, 1.0;
        }

        CorrelationMatrix(const Eigen::VectorXd& coeffs) {

            size_t c = coeffs.size();
            size_t p = static_cast<int>(1.0 + std::sqrt(1.0 + 8.0 * coeffs.size())) / 2.0;

            matrix = Eigen::MatrixXd::Identity(p, p);

            size_t idx = 0;
            for (size_t i = 0; i < p; ++i) {
                for (size_t j = i + 1; j < p; ++j) {
                    double c = coeffs[idx];
                    c = std::abs(c) <= 1.0 ? c : std::numeric_limits<double>::quiet_NaN();
                    matrix(i, j) = c;
                    matrix(j, i) = c; 
                    ++idx;
                }
            }

        }

        CorrelationMatrix(const Eigen::MatrixXd& matrix_): matrix(matrix_){}

        double correlationCoefficient(int i, int j) const { return matrix(i,j); }

        Eigen::MatrixXd eigenvectors() const {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);return solver.eigenvectors(); }

        Eigen::VectorXd eigenvalues() const {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);return solver.eigenvalues(); }

        Eigen::MatrixXd choleskyDecomposition() const {Eigen::LLT<Eigen::MatrixXd> llt(matrix);return llt.matrixL();}

        bool isValid() {

            if (matrix.rows() != matrix.cols()) return false;
            if (!((matrix - matrix.transpose()).norm() < 1e-10)) return false;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
            if (!((solver.eigenvalues().array() > 0).all())) return false;
            for (double v: matrix.diagonal()) if (std::abs(v-1.0)>1e-10) return false;
            for (size_t i = 0; i < matrix.rows(); ++i) {
                for (size_t j = i + 1; j < matrix.rows(); ++j) {
                    if (std::abs(matrix(i, j))-1.0 > 1e-20) return false;
                }
            }
            return true;
        }

        static bool isNumberCoefficientsValid(int n) {
            int p = static_cast<int>(1.0 + std::sqrt(1.0 + 8.0 * n)) / 2.0;
            return !(p * (p - 1) / 2 != n);
        }

    }; 

    struct CovarianceMatrix {
 
        Eigen::MatrixXd matrix;

        CovarianceMatrix(const Eigen::MatrixXd& matrix_): matrix(matrix_){}

        CovarianceMatrix(const CorrelationMatrix& corrMatrix, const Eigen::VectorXd& variances) {

            int n = corrMatrix.matrix.rows();
            matrix = variances.asDiagonal();
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    double cov = corrMatrix.correlationCoefficient(i,j)*std::sqrt(variances(i))*std::sqrt(variances(j));
                    matrix(i, j) = cov; 
                    matrix(j, i) = cov; 
                }
            }
        }

        CorrelationMatrix correlationMatrix() const { 

            size_t n = matrix.rows();
            Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Identity(n,n); 
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    double corrCoeff = correlationCoefficient(i,j);
                    corrMatrix(i, j) = corrCoeff; 
                    corrMatrix(j, i) = corrCoeff; 
                }
            }
            return CorrelationMatrix{corrMatrix};
        }

        double correlationCoefficient(int i, int j) const {return matrix(i,j) / ( std::sqrt(matrix(i,i))*std::sqrt(matrix(j,j)));}

        double variance(int i) const {return matrix(i,i); }

        bool isValid() const {

            size_t n = matrix.rows();
            Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Identity(n,n); 
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    double corrCoeff = correlationCoefficient(i,j);
                    if (std::isnan(corrCoeff)) return false;
                    corrMatrix(i, j) = corrCoeff; 
                    corrMatrix(j, i) = corrCoeff; 
                }
            }
            return CorrelationMatrix{corrMatrix}.isValid();
        }

    };

}