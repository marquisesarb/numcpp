#pragma once 
#include <functional>
#include <numcpp/stats/distributions/normal.hpp>
#include <Eigen/Dense>

namespace numcpp::cfin {

    class BlackScholesEulerSimulation {

        public:
            double T; 
            double dt; 
            double x0;
            Eigen::MatrixXd wienerMatrix; 
            Eigen::MatrixXd localVolatilityMatrix; 
            Eigen::VectorXd timeVector; 
            Eigen::VectorXd driftVector;
            Eigen::MatrixXd xMatrix;

            BlackScholesEulerSimulation(
                double x0_, 
                double T_, 
                const std::function<double(double)>& driftFunction, 
                const std::function<double(double,double)>& localVolatilityFunction, 
                const Eigen::MatrixXd& wienerMatrix_) {

                wienerMatrix = wienerMatrix_;  
                x0 = x0_, T = T_;
                double t=0.0; 
                int N = wienerMatrix_.cols(), M = wienerMatrix_.rows();
                dt = T_/double(N); 
                xMatrix = Eigen::MatrixXd::Zero(M, N);
                localVolatilityMatrix = Eigen::MatrixXd::Zero(M, N);
                timeVector = Eigen::VectorXd::Zero(N);
                driftVector = Eigen::VectorXd::Zero(N);
                double totalDrift = 1.0;
                for (size_t i = 0; i<N; i++) {

                    t += dt; 
                    timeVector(i) = t; 
                    double mut = driftFunction(t); 
                    driftVector(i) = mut;
                    totalDrift *= std::exp(mut*dt);
                    
                    for (size_t z = 0; z<M; z++) {

                        double xt_1 = (i==0) ? x0 : xMatrix(z,i-1);
                        double zt = wienerMatrix(z,i);
                        double sigma = localVolatilityFunction(t,std::log(xt_1/(x0*totalDrift)));
                        double xt = xt_1*std::exp(dt*(mut - .5*sigma*sigma) + sigma*std::sqrt(dt)*zt);
                        xMatrix(z,i) = xt;
                        localVolatilityMatrix(z,i) = sigma;
                    }
                }


            }

            BlackScholesEulerSimulation(
                double x0_,  
                double T_, 
                const std::function<double(double)>& driftFunction, 
                const std::function<double(double,double)>& localVolatilityFunction, 
                int numberPaths, int numberSteps): BlackScholesEulerSimulation(x0_,T_,driftFunction,localVolatilityFunction,generateWienerMatrix(numberPaths,numberSteps)){}

            BlackScholesEulerSimulation(double x0_,  double T_, double sigma, double mu, int numberPaths, int numberSteps): 
            BlackScholesEulerSimulation(x0_,T_,[mu](double t) {return mu;},[sigma](double t,double k) {return sigma;},generateWienerMatrix(numberPaths,numberSteps)){}

            
        private:
            static Eigen::MatrixXd generateWienerMatrix(int numberPaths, int numberSteps) {

                std::random_device rd;
                std::mt19937 g(rd());
                return numcpp::stats::NormalDistribution(0.0,1.0).samples(g,numberPaths,numberSteps);
            }

    }; 

}