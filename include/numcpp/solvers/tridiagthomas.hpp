#pragma once 
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace numcpp::solvers {

    inline Eigen::VectorXd tridiagonalMatrixSolverThomas(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& d) {

        Eigen::VectorXd cp(d.size()); 
        Eigen::VectorXd dp(d.size()); 
        Eigen::VectorXd x(d.size()); 
        cp(0) = A.coeff(0,1)/A.coeff(0,0);
        dp(0) = d(0)/A.coeff(0,0);

        for (size_t i = 1; i<d.size(); i++) {

            
            double a = A.coeff(i,i-1); 
            double b =  A.coeff(i,i); 
            double c = (i==d.size()-1) ? 0 : A.coeff(i,i+1); 

            cp(i) = c/(b - a*cp(i-1)); 
            dp(i) = (d(i)-a*dp(i-1))/(b - a*cp(i-1));
        
        }

        x(d.size()-1) = dp(d.size()-1);

        for (size_t i = d.size()-2; i>0; i--) {
            x(i) = dp(i) - cp(i)*x(i+1);

        }

        return x;

    }

}