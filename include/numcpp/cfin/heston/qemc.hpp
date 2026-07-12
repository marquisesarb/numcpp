#pragma once 
#include <tuple>
#include <numcpp/stats/distributions/normal.hpp>
#include <numcpp/cfin/heston/cf.hpp>

namespace numcpp::cfin {

    struct HestonQuadraticExponentialVarianceSimulation {

        Heston params; 
        double T; 
        double dt; 
        Eigen::MatrixXd uniformMatrix; 
        Eigen::VectorXd timeVector; 
        Eigen::MatrixXd varianceMatrix;
        double phiCritical = 1.5;

        HestonQuadraticExponentialVarianceSimulation(const Heston& params_ ,double T_, const Eigen::MatrixXd& uniformMatrix_) {

            params = params_;
            uniformMatrix = uniformMatrix_; 
            T = T_; 
            dt = T/double(uniformMatrix_.cols()); 
            double t=0.0; 
            int N = uniformMatrix.cols(), M = uniformMatrix.rows();
            varianceMatrix = Eigen::MatrixXd::Zero(M, N);
            timeVector = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i<N; i++) {

                t += dt; 
                timeVector(i) = t; 
                
                for (size_t z = 0; z<M; z++) {

                    double vt_1 = (i==0) ? params.v0 : varianceMatrix(z,i-1);
                    double ut = uniformMatrix(z,i);
                    varianceMatrix(z,i) = computeNextVariance(vt_1, dt, ut);
                }
            }
        }

        HestonQuadraticExponentialVarianceSimulation(const Heston& params,double T_, int numberPaths, int numberSteps): 
        HestonQuadraticExponentialVarianceSimulation(params,T_,generateUniformMatrix(numberPaths, numberSteps)){}

        double computeFirstMoment(double vt, double dt) const {return params.theta + (vt-params.theta)*std::exp(-params.kappa*dt);}

        double computeSecondMoment(double vt, double dt) const {
            double ek = std::exp(-params.kappa*dt); 
            double ftr = params.eta*params.eta * (1-ek) / params.kappa;
            return ftr*(vt*ek + .5*params.theta*(1-ek));
        }

        double computeNextVariance(double vt, double dt, double u) const {

            double m1 = computeFirstMoment(vt, dt);
            double phi = computeSecondMoment(vt,dt)/(m1*m1);

            if (phi <= phiCritical) {

                double two_over_phi = 2.0/phi; 
                double b = std::sqrt(two_over_phi - 1.0 + std::sqrt(two_over_phi*(two_over_phi-1.0)));
                double a = m1/(1.0+b*b);
                double z = numcpp::stats::NormalDistribution(0.0,1.0).invCdf(u);
                return a*(b+z)*(b+z);

            } else {

                double p = (phi-1.0)/(phi+1.0);
                double beta = (1-p)/m1;
                return (u<=p) ? 0.0  : std::log((1-p)/(1-u))/beta;

            }

        }

        static Eigen::MatrixXd generateUniformMatrix(int numberPaths, int numberSteps) {

            std::random_device rd;
            std::mt19937 g(rd());
            return numcpp::stats::UniformDistribution(0.0,1.0).samples(g,numberPaths,numberSteps);
        }

    };

    struct HestonQuadraticExponentialSimulation {

        Heston params; 
        double x0;
        double T;
        double dt;
        Eigen::MatrixXd varUniformMatrix; 
        Eigen::MatrixXd xUniformMatrix; 
        Eigen::VectorXd driftVector;
        Eigen::VectorXd timeVector; 
        Eigen::MatrixXd xMatrix;
        Eigen::MatrixXd varianceMatrix;
        double phiCritical = 1.5;
        double gamma_ = .5; 

        HestonQuadraticExponentialSimulation(
            double x0_, 
            double T_, 
            const Heston& params_,
            const std::function<double(double)>& driftFunction, 
            const Eigen::MatrixXd& uniformMatrix_,
            const Eigen::MatrixXd& uniformMatrix2_) {

            assert(uniformMatrix_.cols() == uniformMatrix2_.cols());
            assert(uniformMatrix_.rows() == uniformMatrix2_.rows());

            params = params_;
            xUniformMatrix = uniformMatrix_; 
            varUniformMatrix = uniformMatrix2_;
            x0 = x0_; 
            T = T_;
            int N = uniformMatrix_.cols();
            int M = uniformMatrix_.rows();
            double t = 0.0;
            double dt = T/double(N);
            xMatrix = Eigen::MatrixXd::Zero(M, N);
            varianceMatrix = Eigen::MatrixXd::Zero(M, N);
            driftVector = Eigen::VectorXd::Zero(N);
            timeVector = Eigen::VectorXd::Zero(N);

            for (size_t i = 0; i<N; i++) {

                t += dt; 
                timeVector(i) = t; 
                double mut = driftFunction(t); 
                driftVector(i) = mut;
                
                for (size_t z = 0; z<M; z++) {

                    double xt_1 = (i==0) ? x0 : xMatrix(z,i-1);
                    double vt_1 = (i==0) ? params.v0 : varianceMatrix(z,i-1);
                    double u1t = xUniformMatrix(z,i); 
                    double u2t = varUniformMatrix(z,i); 
                    double vt = computeNextVariance(vt_1, dt, u2t);
                    double zt = numcpp::stats::NormalDistribution(0.0,1.0).invCdf(u1t);
                    auto [K0,K1,K2,K3,K4] = simulationParams(dt, vt_1, true);
                    double xt = xt_1*std::exp(dt*mut + K0 + K1*vt_1 + K2*vt + std::sqrt(K3*vt_1 + K4*vt)*zt);
                    
                    varianceMatrix(z,i) = vt;
                    xMatrix(z,i) = xt;
                }
            }

        }

        HestonQuadraticExponentialSimulation(
            double x0_, 
            double T_, 
            const Heston& params,
            const std::function<double(double)>& driftFunction, 
            int numberPaths,
            int numberSteps):
        HestonQuadraticExponentialSimulation(x0_,T_,params, driftFunction, generateUniformMatrix(numberPaths, numberSteps), generateUniformMatrix(numberPaths, numberSteps)) {}

        std::tuple<double,double,double,double,double> simulationParams(double dt, double vt, bool martingaleCorrection) const {

            double K0 = -dt*params.kappa*params.rho*params.theta/params.eta; 
            double K1 = (params.kappa*params.rho/params.eta - .5)*gamma_*dt - params.rho/params.eta; 
            double K2 = (params.kappa*params.rho/params.eta - .5)*(1.0-gamma_)*dt + params.rho/params.eta; 
            double K3 = (1.0 - params.rho*params.rho)*gamma_*dt; 
            double K4 = (1.0 - params.rho*params.rho)*(1.0-gamma_)*dt; 

            if (martingaleCorrection) {

                double m1 = computeFirstMoment(vt, dt);
                double phi = computeSecondMoment(vt,dt)/(m1*m1);
                double M;
                double A = K2 + .5*K4;
                if (phi <= phiCritical) {

                    double two_over_phi = 2.0/phi; 
                    double b = std::sqrt(two_over_phi - 1.0 + std::sqrt(two_over_phi*(two_over_phi-1.0)));
                    double a = m1/(1.0+b*b);
                    M = std::exp((A*b*b*a)/(1.0-2.0*A*a))/std::sqrt(1.0-2.0*A*a); 


                } else {

                    double p = (phi-1.0)/(phi+1.0);
                    double beta = (1-p)/m1;
                    M = p+beta*(1-p)/(beta-A);

                }

                K0 = -std::log(M)-(K1+.5*K3)*vt;

            } 

            return {K0,K1,K2,K3,K4};
        }

        double computeFirstMoment(double vt, double dt) const {return params.theta + (vt-params.theta)*std::exp(-params.kappa*dt);}

        double computeSecondMoment(double vt, double dt) const {
            double ek = std::exp(-params.kappa*dt); 
            double ftr = params.eta*params.eta * (1-ek) / params.kappa;
            return ftr*(vt*ek + .5*params.theta*(1-ek));
        }

        double computeNextVariance(double vt, double dt, double u) const {

            double m1 = computeFirstMoment(vt, dt);
            double phi = computeSecondMoment(vt,dt)/(m1*m1);

            if (phi <= phiCritical) {

                double two_over_phi = 2.0/phi; 
                double b = std::sqrt(two_over_phi - 1.0 + std::sqrt(two_over_phi*(two_over_phi-1.0)));
                double a = m1/(1.0+b*b);
                double z = numcpp::stats::NormalDistribution(0.0,1.0).invCdf(u);
                return a*(b+z)*(b+z);

            } else {

                double p = (phi-1.0)/(phi+1.0);
                double beta = (1-p)/m1;
                return (u<=p) ? 0.0  : std::log((1-p)/(1-u))/beta;

            }

        }
        
        static Eigen::MatrixXd generateUniformMatrix(int numberPaths, int numberSteps) {

            std::random_device rd;
            std::mt19937 g(rd());
            return numcpp::stats::UniformDistribution(0.0,1.0).samples(g,numberPaths,numberSteps);
        }


    };
}