#pragma once
#include <functional>
#include <cmath>
#include <numcpp/optim/newtonraphson.hpp>
#include <numcpp/cfin/blackscholes.hpp>


namespace numcpp::cfin {

    namespace bawtools {

        double getDiscountFactor(double t, double r);
        std::function<double(double)> _getTarget(double mu, double K, double t, double sigma, double r, bool isCall);
        std::function<double(double)> _getTargetFirstDerivative(double mu, double KK, double t, double sigma, double r, bool isCall);
        double getQ(double mu, double sigma, double r, double t, bool isCall); 
        double getQInf(double mu, double sigma, double r, bool isCall); 
        std::function<double(double)> _getTargetImpliedVolatiltiy(double p, double F, double K, double t, double r, bool isCall);
        std::function<double(double)> _getTargetFirstDerivativeImpliedVolatiltiy(double p, double F, double K, double t, double r, bool isCall);
        double _getExercisePremium(double F, double mu, double K, double t, double sigma, double r, bool isCall, double E);
        double _getPrice(double F, double mu, double K, double t, double sigma, double r, bool isCall, double E);
        double _getVega(double F, double mu, double K, double t, double sigma, double r, bool isCall, double E1, double E2, double eps);
        double _getImpliedVolatilityNewtonRaphson(double p, double F, double K, double t, double r, bool isCall);
        double europeanPrice(double F, double K, double t, double sigma, double r, bool isCall);
        double europeanDelta(double F, double mu, double K, double t, double sigma, double r, bool isCall);
        double europeanGamma(double F, double mu, double K, double t, double sigma, double r, bool isCall);
        double europeanVega(double F, double K, double t, double sigma, double r, bool isCall);
    }

    inline double baroneAdesiWhaleyOptimalExercisePriceInitialGuess(double mu, double K, double t, double sigma, double r, bool isCall) {

        double putCallFlag = isCall ? 1.0 : -1.0; 
        double qq_inf_ = bawtools::getQInf(mu,sigma,r,isCall);
        double Sinf = K/(1 - 1/qq_inf_);
        double payoff = putCallFlag*(Sinf - K); 
        double h = -putCallFlag*(mu*t+putCallFlag*2*sigma*std::sqrt(t))*K/payoff; 
        return Sinf-putCallFlag*std::exp(h)*payoff;
    }

    inline double baroneAdesiWhaleyOptimalExercisePrice(double mu, double K, double t, double sigma, double r, bool isCall) {

        if (mu >= r and isCall){return 10e9;}
        numcpp::optim::NewtonRaphson res(
            baroneAdesiWhaleyOptimalExercisePriceInitialGuess(mu, K,t,sigma,r,isCall),
            bawtools::_getTarget(mu, K, t, sigma, r, isCall),
            bawtools::_getTargetFirstDerivative(mu, K,t,sigma,r,isCall),
            10, 1e-12
        );
        return res.x;

    }

    inline double baroneAdesiWhaleyExercisePremium(double F, double mu, double K, double t, double sigma, double r, bool isCall) {

        return bawtools::_getExercisePremium(
            F
            ,mu
            ,K
            ,t
            ,sigma
            ,r
            ,isCall
            ,baroneAdesiWhaleyOptimalExercisePrice(mu,K,t,sigma,r,isCall));
    }

    inline double baroneAdesiWhaleyAmericanPrice(double F, double mu, double K, double t, double sigma, double r, bool isCall) {

        return bawtools::_getPrice(
            F
            ,mu
            ,K
            ,t
            ,sigma
            ,r
            ,isCall
            ,baroneAdesiWhaleyOptimalExercisePrice(mu,K,t,sigma,r,isCall));
    }

    inline double baroneAdesiWhaleyImpliedVolatility(double p, double F, double K, double t, double r, bool isCall) {

        if (0.0 >= r and isCall) { 

            return numcpp::cfin::blackScholesImpliedVolatility(
                p/bawtools::getDiscountFactor(t,r)/std::sqrt(F*K)
                , std::log(F/K)
                , t
                , isCall);

        } else {

            return bawtools::_getImpliedVolatilityNewtonRaphson(p,F,K,t,r,isCall);
        }
}



}