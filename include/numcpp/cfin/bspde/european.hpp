#pragma once 
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace numcpp::cfin {


    struct BlackScholesImplicitSchemeEuropeanEngine {


        double dt, dx, smax, S, sigma;
        Eigen::VectorXd terminalSpotPrice; 
        Eigen::MatrixXd optionValues; 

        BlackScholesImplicitSchemeEuropeanEngine(double S_, double K, double mu, double r, double sigma, double T, bool isCall, int N, int M) {

            dt = T/N;
            S = S_;
            terminalSpotPrice = Eigen::VectorXd::Zero(2*M+1);
            optionValues = Eigen::MatrixXd::Zero(2*M+1, N);
            smax = std::max(5.0*std::abs((mu-.5*sigma*sigma))*std::sqrt(dt), std::sqrt(2)*sigma);
            dx = smax * std::sqrt(dt);
            double putCallFlag = isCall ? 1.0 : -1.0;
            terminalSpotPrice = computeTerminalSpotPriceVector(S,dx,M);
            
            
            Eigen::SparseMatrix<double> transitionProbabilityMatrix = computeProbabilityTransitionMatrix(dx, dt, sigma, mu, r, M, smax);

            for (size_t i = 0; i<2*M+1; i++) {
                optionValues(i, N-1) = std::max(putCallFlag*(terminalSpotPrice(i)-K), 0.0);
            }

            for (size_t i = N-1; i>0; i--) {

                optionValues.col(i-1) = transitionProbabilityMatrix*optionValues.col(i);
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

        static Eigen::SparseMatrix<double> computeProbabilityTransitionMatrix(double dx, double dt, double sigma, double mu, double r, int M, double smax) {

            Eigen::SparseMatrix<double> transitionProbabilityMatrix(2*M+1, 2*M+1);
            double p = .5*sigma*sigma/(smax*smax); 
            double d = .5*(mu-.5*sigma*sigma)*std::sqrt(dt)/smax;
            double df = 1.0/(1+r*dt);
            double pu = df*(p+d), pd = df*(p-d), pm = df*(1.0-2.0*p);
            std::vector<Eigen::Triplet<double>> coefficients;

            coefficients.push_back({Eigen::Triplet<double>(0,0,pu)});
            coefficients.push_back({Eigen::Triplet<double>(0,1,2.0*pm-pu)});
            coefficients.push_back({Eigen::Triplet<double>(0,2,2.0*pd-pm)});
            coefficients.push_back({Eigen::Triplet<double>(0,3,pd)});
            coefficients.push_back({Eigen::Triplet<double>(2*M,2*M,df)});

            for (size_t i = 1; i<2*M; i++) {

                coefficients.push_back({Eigen::Triplet<double>(i,i-1,pu)});
                coefficients.push_back({Eigen::Triplet<double>(i,i,pm)});
                coefficients.push_back({Eigen::Triplet<double>(i,i+1,pd)});
            }

            transitionProbabilityMatrix.setFromTriplets(coefficients.begin(), coefficients.end());
            return transitionProbabilityMatrix;
        }

        double price() const {

            int M = (terminalSpotPrice.size()-1)/2;
            return optionValues.col(0)(M);
        }

    };


    struct BlackScholesImplicitSchemeEuropeanEngine2 {


        double dt, dx, smax, S, sigma;
        Eigen::VectorXd terminalSpotPrice; 
        Eigen::MatrixXd optionValues; 

        BlackScholesImplicitSchemeEuropeanEngine2(double S_, double K, double mu, double r, double sigma, double T, bool isCall, int N, int M) {

            dt = T/N;
            S = S_;
            terminalSpotPrice = Eigen::VectorXd::Zero(2*M+1);
            optionValues = Eigen::MatrixXd::Zero(2*M+1, N);
            smax = std::max(5.0*std::abs((mu-.5*sigma*sigma))*std::sqrt(dt), std::sqrt(2)*sigma);
            dx = smax * std::sqrt(dt);
            double putCallFlag = isCall ? 1.0 : -1.0;
            terminalSpotPrice = computeTerminalSpotPriceVector(S,dx,M);
            Eigen::SparseMatrix<double> transitionProbabilityMatrix = computeProbabilityTransitionMatrix(dx, dt, sigma, mu, r, M, smax);

            for (size_t i = 0; i<2*M+1; i++) {
                optionValues(i, N-1) = std::max(putCallFlag*(terminalSpotPrice(i)-K), 0.0);
            }

            for (size_t i = N-1; i>0; i--) {

                optionValues.col(i-1) = transitionProbabilityMatrix*optionValues.col(i);
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

        static Eigen::SparseMatrix<double> computeProbabilityTransitionMatrix(double dx, double dt, double sigma, double mu, double r, int M, double smax) {

            Eigen::SparseMatrix<double> transitionProbabilityMatrix(2*M+1, 2*M+1);
            double p = .5*sigma*sigma/(smax*smax); 
            double d = .5*(mu-.5*sigma*sigma)*std::sqrt(dt)/smax;
            double df = 1.0/(1+r*dt);
            double pu = df*(p+d), pd = df*(p-d), pm = df*(1.0-2.0*p);
            std::vector<Eigen::Triplet<double>> coefficients;

            coefficients.push_back({Eigen::Triplet<double>(0,0,pu)});
            coefficients.push_back({Eigen::Triplet<double>(0,1,2.0*pm-pu)});
            coefficients.push_back({Eigen::Triplet<double>(0,2,2.0*pd-pm)});
            coefficients.push_back({Eigen::Triplet<double>(0,3,pd)});
            coefficients.push_back({Eigen::Triplet<double>(2*M,2*M,df)});

            for (size_t i = 1; i<2*M; i++) {

                coefficients.push_back({Eigen::Triplet<double>(i,i-1,pu)});
                coefficients.push_back({Eigen::Triplet<double>(i,i,pm)});
                coefficients.push_back({Eigen::Triplet<double>(i,i+1,pd)});
            }

            transitionProbabilityMatrix.setFromTriplets(coefficients.begin(), coefficients.end());
            return transitionProbabilityMatrix;
        }

        double price() const {

            int M = (terminalSpotPrice.size()-1)/2;
            return optionValues.col(0)(M);
        }

    };

}