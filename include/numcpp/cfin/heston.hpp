#pragma once 
#include <complex>
#include <numcpp/cfin/lewis.hpp>

namespace numcpp::cfin {


    struct Heston {

        double v0, theta, rho, eta, kappa;

        std::function<std::complex<double>(std::complex<double>)> characteristicFunction(double t) const {

            auto cf = [*this,t](std::complex<double> u){ 
                std::complex<double> i(0.0, 1.0);
                std::complex<double> beta = kappa - rho*eta*i*u;
                std::complex<double> d = std::sqrt(beta*beta + eta*eta*(i*u + u*u));
                std::complex<double> g = (beta - d)/(beta + d); 
                std::complex<double> D = (beta - d)*(1.0-std::exp(-d*t))/(eta*eta*(1.0-g*std::exp(-d*t)));
                std::complex<double> C = kappa*(t*(beta-d) - 2.0*std::log((1.0-g*std::exp(-d*t))/(1.0-g)))/(eta*eta); 
                return std::exp(C*theta+D*v0); 
            };

            return cf;
        }

        double europeanNormalizedPriceGaussLaguerre(double x, double t, bool isCall, const numcpp::gquad::GaussLaguerreQuadrature& gaussLaguerreQuadratures) {

            return lewisEuropeanNormalizedPriceGaussLaguerre(x,isCall,characteristicFunction(t),gaussLaguerreQuadratures);
        }

        double impliedVolatilityGaussLaguerre(double x, double t, bool isCall, const numcpp::gquad::GaussLaguerreQuadrature& gaussLaguerreQuadratures) {

            return lewisImpliedVolatilityGaussLaguerre(x, t, isCall, characteristicFunction(t),gaussLaguerreQuadratures);

        }

        double impliedDensityGaussLaguerre(double x, double t, const numcpp::gquad::GaussLaguerreQuadrature& gaussLaguerreQuadratures) {

            return lewisImpliedDensityGaussLaguerre(x, characteristicFunction(t), gaussLaguerreQuadratures);
        }


    };
}