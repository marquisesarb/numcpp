#pragma once 
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <numcpp/solvers/tridiagthomas.hpp>

namespace numcpp::cfin {


    struct BlackScholesImplicitSchemeEuropeanEngine {

        Eigen::VectorXd terminalSpotPrice; 
        Eigen::MatrixXd optionValues; 

        BlackScholesImplicitSchemeEuropeanEngine(double S_, double K, double mu, double r, double sigma, double T, bool isCall, int N, int M) {

            double dt = T/N;
            double S = S_;
            terminalSpotPrice = Eigen::VectorXd::Zero(2*M+1);
            optionValues = Eigen::MatrixXd::Zero(2*M+1, N);
            double xmax = (mu - .5*sigma*sigma)*T + 5.0*std::sqrt(T)*sigma;  
            double dx = xmax/double(M);
            double putCallFlag = isCall ? 1.0 : -1.0;
            terminalSpotPrice = computeTerminalSpotPriceVector(S,dx,M);

            Eigen::SparseMatrix<double> transitionMatrix = computeMatrix(dx, dt, sigma, mu, r, M);

            for (size_t i = 0; i<2*M+1; i++) {

                optionValues(i, N-1) = std::max(putCallFlag*(terminalSpotPrice(i)-K), 0.0);
            }

            for (size_t i = N-1; i>0; i--) {

                optionValues.col(i-1).segment(1,2*M) = solvers::tridiagonalMatrixSolverThomas(transitionMatrix, optionValues.col(i).segment(1,2*M)); 
                optionValues.col(i-1)(0) = 2*optionValues.col(i-1)(1) - optionValues.col(i-1)(2);
            }

        }

        static Eigen::VectorXd computeTerminalSpotPriceVector(double S, double dx, int M) {

            Eigen::VectorXd logMoneynessVector_ = Eigen::VectorXd::Zero(2*M+1);
            Eigen::VectorXd terminalSpotPrice_ = Eigen::VectorXd::Zero(2*M+1);
            size_t mid = M;
            logMoneynessVector_(mid) = 0.0;
            terminalSpotPrice_(mid) = S;
            for (size_t i = 1; i<=mid; i++) {
                logMoneynessVector_(mid-i) = logMoneynessVector_(mid-i+1)+dx;
                logMoneynessVector_(mid+i) = logMoneynessVector_(mid+i-1)-dx;
                terminalSpotPrice_(mid-i) = S*std::exp(logMoneynessVector_(mid-i));
                terminalSpotPrice_(mid+i) = S*std::exp(logMoneynessVector_(mid+i));

            }
            return terminalSpotPrice_;
        }

        static Eigen::SparseMatrix<double> computeMatrix(double dx, double dt, double sigma, double mu, double r, int M) {

            Eigen::SparseMatrix<double> transitionProbabilityMatrix(2*M, 2*M);
            std::vector<Eigen::Triplet<double>> coefficients;
            double mid = 1.0+dt*r + sigma*sigma*dt/(dx*dx); 
            double up = -.5*sigma*sigma*dt/(dx*dx) - .5*(mu - .5*sigma*sigma)*dt/dx; 
            double down = -.5*sigma*sigma*dt/(dx*dx) + .5*(mu - .5*sigma*sigma)*dt/dx; 

            coefficients.push_back({Eigen::Triplet<double>(0,0,mid)});
            coefficients.push_back({Eigen::Triplet<double>(0,1,down)});
            coefficients.push_back({Eigen::Triplet<double>(2*M-1,2*M-1,1.0+dt*r)});

            for (size_t i = 1; i<2*M-1; i++) {

                coefficients.push_back({Eigen::Triplet<double>(i,i-1,up)});
                coefficients.push_back({Eigen::Triplet<double>(i,i,mid)});
                coefficients.push_back({Eigen::Triplet<double>(i,i+1,down)});
            }

            transitionProbabilityMatrix.setFromTriplets(coefficients.begin(), coefficients.end());
            return transitionProbabilityMatrix;
        }

        double price() const {

            return optionValues.col(0)(size_t((terminalSpotPrice.size()-1)/2));
        }

    };

}