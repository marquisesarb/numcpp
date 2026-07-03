#pragma once 
#include <functional>
#include <limits>
#include <utility>

namespace numcpp::optim {

    struct NewtonRaphson {

        double x, fx, dfx, toleranceThreshold, zeroThreshold;
        int iter,  maxIter;
        bool converged;
 
        NewtonRaphson(
            double x0
            , const std::function<std::pair<double,double>(double)>& fdf
            , int maxIter_ = 100
            , double toleranceThreshold_ = 1e-12
            , double zeroThreshold_ = 1e-20) {

            zeroThreshold = zeroThreshold_; 
            toleranceThreshold = toleranceThreshold_; 
            maxIter = maxIter_;
            x = x0;
            double xStep = 0;
            iter = 0;
            while (true) {

                std::pair<double,double> fwDerivAtX = fdf(x);
                fx = fwDerivAtX.first;
                dfx = fwDerivAtX.second;
                iter += 1;
                if (std::abs(fx) < toleranceThreshold){converged = true; break;}
                
                if (std::abs(dfx) < zeroThreshold) {

                    x=std::numeric_limits<double>::quiet_NaN();
                    fx=std::numeric_limits<double>::quiet_NaN();
                    dfx=std::numeric_limits<double>::quiet_NaN();
                    converged = false;
                    break;
                }
                
                double xNew = x - fx / dfx;
                xStep = std::abs(xNew - x);
                x = xNew;

                if (xStep < toleranceThreshold) {converged = true; break;}
                if (iter==maxIter) {converged = false; break;}
                
                
            }

        }

        NewtonRaphson(
            double x0
            , const std::function<double(double)>& f
            , const std::function<double(double)>& df
            , int maxIter_ = 100
            , double toleranceThreshold_ = 1e-12
            , double zeroThreshold_ = 1e-20): 
        NewtonRaphson(x0, [f,df](double x){return std::make_pair(f(x), df(x));}, maxIter_,toleranceThreshold_,zeroThreshold_){}



    };
}