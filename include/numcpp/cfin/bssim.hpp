#pragma once 
#include <functional>
#include <numcpp/stats/distributions/normal.hpp>
#include <Eigen/Dense>

namespace numcpp::cfin {

    struct BlackScholesSimulation {

        double sigma; 
        double T; 
        double dt; 
        double x0;
        Eigen::MatrixXd wienerMatrix; 
        Eigen::VectorXd timeVector; 
        Eigen::VectorXd driftVector;
        Eigen::MatrixXd xMatrix;

        BlackScholesSimulation(
            double x0_, 
            double sigma_, 
            double T_, 
            const std::function<double(double)>& driftFunction, 
            const Eigen::MatrixXd& wienerMatrix_) {

            wienerMatrix = wienerMatrix_;  
            sigma = sigma_, x0 = x0_, T = T_;
            double t=0.0; 
            int N = wienerMatrix_.cols(), M = wienerMatrix_.rows();
            dt = T_/double(N); 
            xMatrix = Eigen::MatrixXd::Zero(M, N);
            timeVector = Eigen::VectorXd::Zero(N);
            driftVector = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i<N; i++) {

                t += dt; 
                timeVector(i) = t; 
                double mut = driftFunction(t); 
                driftVector(i) = mut;
                
                for (size_t z = 0; z<M; z++) {

                    double xt_1 = (i==0) ? x0 : xMatrix(z,i-1);
                    double zt = wienerMatrix(z,i);
                    double xt = xt_1*std::exp(dt*(mut - .5*sigma*sigma) + sigma*std::sqrt(dt)*zt);
                    xMatrix(z,i) = xt;
                }
            }


        }

        BlackScholesSimulation(
            double x0_, 
            double sigma_, 
            double T_, 
            const std::function<double(double)>& driftFunction, 
            int numberPaths, int numberSteps): BlackScholesSimulation(x0_,sigma_,T_,driftFunction,generateWienerMatrix(numberPaths,numberSteps)){}

        static Eigen::MatrixXd generateWienerMatrix(int numberPaths, int numberSteps) {

            std::random_device rd;
            std::mt19937 g(rd());
            return numcpp::stats::NormalDistribution(0.0,1.0).samples(g,numberPaths,numberSteps);
        }

    }; 

}