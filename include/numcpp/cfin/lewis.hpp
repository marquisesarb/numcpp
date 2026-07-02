#pragma once
#include <complex>
#include <numcpp/gquad/gausslaguerre.hpp> 
#include <numcpp/cfin/blackscholes.hpp>

namespace numcpp::cfin {

    inline double lewisEuropeanNormalizedPriceGaussLaguerre(
        double x
        , bool isCall
        , const std::function<std::complex<double>(std::complex<double>)>& cf
        , const numcpp::gquad::GaussLaguerreQuadrature& gaussLaguerreQuadratures) {

        constexpr double LAGUERRE_WEIGHT_MIN = 1e-20;
        constexpr double INTEGRAND_REAL_PART_MIN = 1e-20;

        size_t n = gaussLaguerreQuadratures.size();
        double sum = 0.0, u, w;
        std::complex<double> integrand;
        std::complex<double> i(0.0, 1.0);
        for (size_t j = 0; j < n; ++j) {
            
            u = gaussLaguerreQuadratures.roots[j]; w = gaussLaguerreQuadratures.weights[j];
            integrand = std::exp(i * u * x) * cf(std::complex<double>(u, -0.5));
            if (abs(integrand.real()) < INTEGRAND_REAL_PART_MIN or abs(w)< LAGUERRE_WEIGHT_MIN) {break;}
            else {sum += w * std::exp(u) * integrand.real() / (u * u + 0.25);}
            
        }
        double call = std::exp(x / 2.0) - sum / numcpp::cfin::letsberational::PI;
        return isCall ? call : call - numcpp::cfin::letsberational::_getNormalizedIntrisicValue(x, true) + numcpp::cfin::letsberational::_getNormalizedIntrisicValue(x, false);
    }

    inline double lewisImpliedVolatilityGaussLaguerre(
        double x
        , double t
        , bool isCall
        , const std::function<std::complex<double>(std::complex<double>)>& cf
        , const numcpp::gquad::GaussLaguerreQuadrature& gaussLaguerreQuadratures) {

        double price = lewisEuropeanNormalizedPriceGaussLaguerre(x, isCall,cf, gaussLaguerreQuadratures);
        return numcpp::cfin::blackScholesImpliedVolatility(price, x, t, isCall);
    } 

    inline double lewisImpliedDensityGaussLaguerre(
        double x
        , const std::function<std::complex<double>(std::complex<double>)>& cf
        , const numcpp::gquad::GaussLaguerreQuadrature& gaussLaguerreQuadratures) {

        constexpr double LAGUERRE_WEIGHT_MIN = 1e-20;
        constexpr double INTEGRAND_REAL_PART_MIN = 1e-20;
        
        size_t n = gaussLaguerreQuadratures.size();
        double sum = 0.0, u, w;
        std::complex<double> integrand;
        std::complex<double> i(0.0, 1.0);
        for (size_t j = 0; j < n; ++j) {
            
            u = gaussLaguerreQuadratures.roots[j]; w = gaussLaguerreQuadratures.weights[j];
            integrand = std::exp(-i * u * x) * cf(std::complex<double>(u, 0.0));
            if (abs(integrand.real()) < INTEGRAND_REAL_PART_MIN or abs(w)< LAGUERRE_WEIGHT_MIN) {break;}
            else {sum += w * std::exp(u) * integrand.real();}
        }
        return sum / numcpp::cfin::letsberational::PI;
    }

}