#pragma once 
#include <functional>
#include <Eigen/Dense>

namespace numcpp::gquad {


    struct GaussLaguerreQuadrature {

        Eigen::VectorXd roots, weights;

        GaussLaguerreQuadrature(int n) {

            int points = std::max(n,5);
            weights = Eigen::VectorXd::Zero(points);
            roots = Eigen::VectorXd::Zero(points);
            Eigen::MatrixXd J = Eigen::MatrixXd::Zero(points, points);

            for (int k = 0; k < points; k++)
            {
                J(k, k) = 2.0 * k + 1.0;
            }

            for (int k = 0; k < points - 1; k++)
            {
                double b = k + 1.0;
                J(k, k + 1) = b;
                J(k + 1, k) = b;
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(J);
            roots = solver.eigenvalues();
            Eigen::MatrixXd V = solver.eigenvectors();
            for (int i = 0; i < points; i++)
            {
                weights[i] = std::pow(V(0, i), 2);
            }
        }

        double integrate(const std::function<double(double)>& f) const {

            double integral = 0.0;
            for (size_t i = 0; i < roots.size(); i++){integral += weights[i] * f(roots[i]);}
            return integral;
        }

        int size() const {return roots.size();}
    };

}