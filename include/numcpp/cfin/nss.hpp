#pragma once
#include <cmath>
#include <Eigen/Dense>

namespace numcpp::cfin {

    struct NelsonSiegel {

        double beta0, beta1, beta2, tau; 

        double spotRate(double t) const {

            double f1 = tau*(1-std::exp(-t/tau))/t; 
            double f2 = tau*(1-std::exp(-t/tau))/t - std::exp(-t/tau);
            return beta0 + beta1 * f1 + beta2 * f2;
        }

        double instantaneousForwardRate(double t) const {

            double f1 = std::exp(-t/tau);
            double f2 = t*std::exp(-t/tau)/tau;
            return beta0 + beta1  * f1 + beta2 * f2;
        }

        double derivativeInstantaneousForwardRate(double t) const {

            return std::exp(-t/tau)*(beta2*(1-t/tau) - beta1);
        }

        double instantaneousSpotRate() const {return beta0 + beta1;}

        Eigen::VectorXd spotRateGradients(double t) {

            Eigen::VectorXd gradient(4); 
            gradient(0) = 1.0; 
            gradient(1) = tau*(1-std::exp(-t/tau))/t; 
            gradient(2) = tau*(1-std::exp(-t/tau))/t - std::exp(-t/tau);
            double dtf = (1.0-std::exp(-t/tau))/t - std::exp(-t/tau)/tau;
            gradient(3) = beta1*dtf + beta2*(dtf - t*std::exp(-t/tau)/(tau*tau));
            return gradient;
        }

        Eigen::VectorXd instantaneousForwardRateGradients(double t) {

            Eigen::VectorXd gradient(4); 
            gradient(0) = 1.0; 
            gradient(1) = std::exp(-t/tau); 
            gradient(2) = t*std::exp(-t/tau)/tau;         
            gradient(3) = t*std::exp(-t/tau)*(beta1+beta2*(t/tau - 1.0))/(tau*tau);
            return gradient;
        }
    };

    struct Svensson {

        double beta0, beta1, beta2, beta3, tau1, tau2; 

        double spotRate(double t) const {

            double f1 = tau1*(1-std::exp(-t/tau1))/t; 
            double f2 = tau1*(1-std::exp(-t/tau1))/t - std::exp(-t/tau1);
            double f3 = tau2*(1-std::exp(-t/tau2))/t - std::exp(-t/tau2);
            return beta0 + beta1 * f1 + beta2 * f2 + beta3*f3;
        }

        double instantaneousForwardRate(double t) const {

            double f1 = std::exp(-t/tau1);
            double f2 = t*std::exp(-t/tau1)/tau1;
            double f3 = t*std::exp(-t/tau2)/tau2;
            return beta0 + beta1  * f1 + beta2 * f2 + beta3*f3;
        }

        double derivativeInstantaneousForwardRate(double t) const {

            return std::exp(-t/tau1)*(beta2*(1-t/tau1) - beta1) + beta3*std::exp(-t/tau2)*(1-t/tau2)/tau2;
        }

        double instantaneousSpotRate() const {return beta0 + beta1;}

        Eigen::VectorXd spotRateGradients(double t) const {

            Eigen::VectorXd gradient(6); 
            gradient(0) = 1.0; 
            gradient(1) = tau1*(1-std::exp(-t/tau1))/t; 
            gradient(2) = tau1*(1-std::exp(-t/tau1))/t - std::exp(-t/tau1);
            gradient(3) = tau2*(1-std::exp(-t/tau2))/t - std::exp(-t/tau2);
            double dtau1f = (1.0-std::exp(-t/tau1))/t - std::exp(-t/tau1)/tau1;
            gradient(4) = beta1*dtau1f + beta2*(dtau1f - t*std::exp(-t/tau1)/(tau1*tau1));
            double dtau2f = (1.0-std::exp(-t/tau2))/t - std::exp(-t/tau2)/tau2;
            gradient(5) = beta3*(dtau2f - t*std::exp(-t/tau2)/(tau2*tau2));
            return gradient;
        }

        Eigen::VectorXd instantaneousForwardRateGradients(double t) const {

            Eigen::VectorXd gradient(6); 
            gradient(0) = 1.0; 
            gradient(1) = std::exp(-t/tau1); 
            gradient(2) = t*std::exp(-t/tau1)/tau1;     
            gradient(3) = t*std::exp(-t/tau2)/tau2;         
            gradient(4) = t*std::exp(-t/tau1)*(beta1+beta2*(t/tau1 - 1.0))/(tau1*tau1);
            gradient(5) = beta3*t*std::exp(-t/tau2)*(t/(tau2*tau2) - tau2)/(tau2*tau2);
            return gradient;
        }    
    };

}