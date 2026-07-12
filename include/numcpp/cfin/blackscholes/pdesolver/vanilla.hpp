#pragma once 
#include <numcpp/cfin/blackscholes/pdesolver/interface.hpp>
#include <numcpp/solvers/tridiagthomas.hpp>

namespace numcpp::cfin {

    class BlackScholesPDESolverVanilla final: public BlackScholesPDESolverInterface {

        public:
            BlackScholesPDESolverVanilla(double S_, double K, double T, double mu, double r, double sigma, bool isCall, bool isAmerican, int N, int M, bool useCrankNicolson, double numberSigma = 5):
            BlackScholesPDESolverInterface(N,M) {

                double dt = T/N;
                double S = S_;
                double xmax = (mu - .5*sigma*sigma)*T + std::abs(numberSigma)*std::sqrt(T)*sigma;  
                double dx = xmax/double(M);
                double putCallFlag = isCall ? 1.0 : -1.0;

                computeLogMoneynessVector(S,dx);
                computeSpotPriceMatrix(S,dt,[mu](double t) {return mu;});
                computeTerminalPayoff(putCallFlag,K);
                
                if (useCrankNicolson) {

                    std::pair<Eigen::SparseMatrix<double>,Eigen::SparseMatrix<double>> mats = flatSigmaCrankNicolsonMatrixes(dx, dt, sigma, 0.0, r, M);

                    for (size_t i = N-1; i>0; i--) {

                        optionValueMatrix.col(i-1) = solvers::tridiagonalMatrixSolverThomas(mats.second, mats.first*optionValueMatrix.col(i));
                        optionValueMatrix.col(i-1)(0) = 2*optionValueMatrix.col(i-1)(1) - optionValueMatrix.col(i-1)(2); 
                        if (isAmerican) computeEarlyExercise(i-1, K, putCallFlag);
                    }

                } else {

                    Eigen::SparseMatrix<double> implicitMatrix = flatSigmaImplicitMatrix(dx, dt, sigma, 0.0, r, M);
                    
                    for (size_t i = N-1; i>0; i--) {

                        optionValueMatrix.col(i-1) = solvers::tridiagonalMatrixSolverThomas(implicitMatrix, optionValueMatrix.col(i));
                        optionValueMatrix.col(i-1)(0) = 2*optionValueMatrix.col(i-1)(1) - optionValueMatrix.col(i-1)(2); 
                        if (isAmerican) computeEarlyExercise(i-1, K, putCallFlag);
                    }


                }


            }

            BlackScholesPDESolverVanilla(
                double S_, double K, double T, bool isCall, bool isAmerican, 
                const std::function<double(double)>& driftFunction, 
                const std::function<double(double,double)>& localVolatilityFunction, 
                const std::function<double(double)>& discountRateFunction,
                int N, int M, bool useCrankNicolson, double numberSigma = 5):
            BlackScholesPDESolverInterface(N,M) {

                double dt = T/N;
                double S = S_;
                double atmLocalVol = localVolatilityFunction(0.0,T);
                double xmax = - .5*atmLocalVol*atmLocalVol*T + std::abs(numberSigma)*std::sqrt(T)*atmLocalVol;  
                double dx = xmax/double(M);
                double putCallFlag = isCall ? 1.0 : -1.0;
                double t = T;

                computeLogMoneynessVector(S,dx);
                computeTerminalPayoff(putCallFlag,K);
                computeSpotPriceMatrix(S,dt,driftFunction);

                if (useCrankNicolson) {

                    for (size_t i = N-1; i>0; i--) {
                        t -= dt; 
                        std::pair<Eigen::SparseMatrix<double>,Eigen::SparseMatrix<double>> mats = localCrankNicolsonMatrixes(t, dx,dt, M, driftFunction, localVolatilityFunction, discountRateFunction);
                        optionValueMatrix.col(i-1) = solvers::tridiagonalMatrixSolverThomas(mats.second, mats.first*optionValueMatrix.col(i));
                        optionValueMatrix.col(i-1)(0) = 2*optionValueMatrix.col(i-1)(1) - optionValueMatrix.col(i-1)(2); 
                        if (isAmerican) computeEarlyExercise(i-1, K, putCallFlag);
                    }

                } else {
                    
                    for (size_t i = N-1; i>0; i--) {
                        t -= dt; 
                        Eigen::SparseMatrix<double> implicitMatrix = localImplicitMatrix(t, dx,dt, M, driftFunction, localVolatilityFunction, discountRateFunction);
                        optionValueMatrix.col(i-1) = solvers::tridiagonalMatrixSolverThomas(implicitMatrix, optionValueMatrix.col(i));
                        optionValueMatrix.col(i-1)(0) = 2*optionValueMatrix.col(i-1)(1) - optionValueMatrix.col(i-1)(2); 
                        if (isAmerican) computeEarlyExercise(i-1, K, putCallFlag);
                    }


                }


            }

        private:
        
            void computeTerminalPayoff(double putCallFlag, double K) {

                for (size_t i = 0; i<spotPriceMatrix.rows(); i++) {
                    optionValueMatrix(i, optionValueMatrix.cols()-1) = std::max(putCallFlag*(spotPriceMatrix(i,spotPriceMatrix.cols()-1)-K), 0.0);
                }
            }

            void computeEarlyExercise(size_t index, double K, double putCallFlag) {

                for (size_t m = 0; m<spotPriceMatrix.rows(); m++) {

                    double optionValue = optionValueMatrix.col(index)(m);
                    double spotValue = spotPriceMatrix.col(index)(m);
                    optionValueMatrix.col(index)(m) = std::max(optionValue, std::max(putCallFlag*(spotValue-K), 0.0));
                }


            }
    };

}