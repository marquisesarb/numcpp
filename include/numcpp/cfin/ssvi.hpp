#pragma once 
#include <cmath>
#include <functional>
#include <Eigen/Dense>
#include <numcpp/cfin/blackscholes/cf.hpp>

namespace numcpp::cfin {

    struct SSVIPowerLaw {

        double rho, nu, g; 
        std::function<double(double)> atmTotalVarianceTermStructure; 

        SSVIPowerLaw(double rho_, double nu_, double g_, std::function<double(double)> atmTotalVarianceTermStructure_): rho(rho_), nu(nu_), g(g_), atmTotalVarianceTermStructure(atmTotalVarianceTermStructure_){}

        virtual ~SSVIPowerLaw() = default;

        double atmTotalVariance(double t) const {return atmTotalVarianceTermStructure(t);}

        double totalVariance(double k, double t) const {

            double w = atmTotalVariance(t); 
            double pl = nu*std::pow(w,-g); 
            double sqterm = std::sqrt((pl*k +rho)*(pl*k +rho) + (1.0- rho*rho));
            return .5*w*(1.0 + pl*rho*k + sqterm);
        }

        double totalVarianceSkew(double k, double t) const {

            double w = atmTotalVariance(t);  
            double pl = nu*std::pow(w,-g); 
            double sqterm = std::sqrt((pl*k +rho)*(pl*k +rho) + (1.0- rho*rho));
            return .5*w*(rho*pl + pl*(pl*k+rho)/sqterm);
        }

        double totaVarianceSecondDerivativeK(double k, double t) const {

            double w = atmTotalVariance(t);  
            double pl = nu*std::pow(w,-g); 
            double sqterm = std::sqrt((pl*k +rho)*(pl*k +rho) + (1.0- rho*rho)); 
            return .5*pl*pl*w*(1-rho*rho)/std::pow(sqterm, 3.0); 
        }

        double totaVarianceFirstDerivativeT(double k, double t, double eps = 1e-5) const {

            return (totalVariance(k,t+eps)-totalVariance(k,t-eps))/(2.0*eps);
        }

        double atmTotalVarianceSkew(double k, double t) const {

            return totalVarianceSkew(0.0, t);
        }

        double densityFactor(double k, double t) const {

            double w = totalVariance(k, t);
            double wk = totalVarianceSkew(k,t);
            double wkk = totaVarianceSecondDerivativeK(k,t);
            return (1.0-0.5*k*wk/w)*(1.0-0.5*k*wk/w)-0.25*wk*wk*(0.25+1/w) + 0.5*wkk;
        }

        double localVariance(double k, double t, double eps = 1e-5) const {

            return totaVarianceFirstDerivativeT(k,t,eps)/densityFactor(k,t);
        }

        double gradientRho(double k, double t) const {

            return 0.0;
        }

        double gradientNu(double k, double t) const {

            return 0.0;
        }

        double gradientGamma(double k, double t) const {

            return 0.0;
        }

        double normalizedPrice(double k, double t, bool isCall) const {

            return numcpp::cfin::blackScholesEuropeanNormalizedPrice(k, t, std::sqrt(totalVariance(k, t)/t), isCall);
        }


    };


    struct SSVIPowerLawFlatATM : public SSVIPowerLaw {

        SSVIPowerLawFlatATM(double rho_, double nu_, double g_, double flatVol): SSVIPowerLaw(rho_,nu_,g_,[flatVol](double t){return flatVol*flatVol*t;}){}
    };
}