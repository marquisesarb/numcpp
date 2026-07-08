#pragma once 
#include <cmath>

namespace nucmpp::cfin {


    inline double sabrImpliedVolatility(double F, double K, double T, double v, double alpha, double beta, double rho, double shift = 0.0) {


        double g = (beta - 1.0)*(beta - 1.0)*std::pow(F+shift,beta-1.0)*std::pow(K+shift,beta-1.0)*alpha*alpha;
        double c = (std::pow(F+shift,beta-1.0)-std::pow(K+shift,beta-1.0))*v/(alpha*(1.0-beta)); 
        double x = std::log((std::sqrt(1.0-2.0*rho*c+c*c) - rho + c)/(1-rho))/v; 
        double factorT = g + .25*rho*v*alpha*std::pow(F+shift,.5*(beta-1.0))*std::pow(K+shift,.5*(beta-1.0)) + (2.0-3.0*rho*rho)*v*v/24.0;
        return std::log(std::pow(F+shift,beta-1.0)/std::pow(K+shift,beta-1.0))*(1.0+T*factorT);
    }
    

}