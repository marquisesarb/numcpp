#pragma once 
#include <functional>
#include <numcpp/optim/newtonraphson.hpp>

namespace numcpp::cfin {

     namespace letsberational {

        static constexpr double H_LARGE = -10.0;
        static constexpr double T_SMALL = 0.21;
        constexpr double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;   
        constexpr double ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887; 
        constexpr double PI = 3.14159265358979323846;
        constexpr double SQRT_TWO_PI = 2.50662827463100050242;

        double _getNormalizedIntrisicValue(double x, bool isCall);
        double _getCallPriceRegion1(double h, double t); 
        double _getCallPriceRegion2(double h, double t); 
        double _getCallPriceRegion3(double h, double t); 
        double _getCallPriceRegion4(double h, double t); 
        double _getCallPrice(double x, double normalizedSigma);
        double _getVega(double x, double normalizedSigma);
        double _getVolga(double x, double normalizedSigma);
        double _getNormalizedPrice(double x, double normalizedSigma, bool isCall);
        double _getNormalizedVolatility(double normalizedPrice, double x, bool isCall);

        double _getRationalCubicInterpolate(double x, double x0, double x1,double y0, double y1,double dy0, double dy1,double r); 

        double _getRightAsymptoticR(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy0);
        double _getRightAsymptotic(double sigma); 
        double _getFirstDerivativeRightAsymptotic(double x, double sigma);
        double _getSecondDerivativeRightAsymptotic(double x, double sigma);

        double _getLeftAsymptoticR(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy1);
        double _getLeftAsymptotic(double x, double sigma); 
        double _getFirstDerivativeLeftAsymptotic(double x, double sigma); 
        double _getSecondDerivativeLeftAsymptotic(double x, double sigma); 
        double _getLeftAsymptoticZ(double x, double sigma);

        std::tuple<double, double, double, double, double, double, double, double, double, double, bool> _getInitialData(double beta, double x, double isCall);
        double _getInitialGuess(double beta, double x, bool isCall); 
        std::function<double(double)> _getTarget(double beta, double x, double bLower, double bUpper);
        std::function<double(double)> _getTargetFirstDerivative(double beta, double x, double bLower, double bUpper); 
        double _getNormalizedPrice(double x, double normalizedSigma, bool isCall);

    }

    inline double blackScholesEuropeanNormalizedPrice(double x, double t, double sigma, bool isCall) {

        return letsberational::_getNormalizedPrice(x, sigma * std::sqrt(t), isCall);
    }

    inline double blackScholesImpliedVolatility(double normalizedPrice, double x, double t,  bool isCall) {

        return letsberational::_getNormalizedVolatility(normalizedPrice, x, isCall)/std::sqrt(t);
    }


}