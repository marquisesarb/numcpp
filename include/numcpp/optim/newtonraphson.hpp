#pragma once 
#include <functional>
#include <limits>
#include <utility>

namespace numcpp::optim {

    struct NewtonRaphson {

        double x, fx, dfx, toleranceThreshold, zeroThreshold;
        int iter,  maxIter;
        
        NewtonRaphson(
            double x0
            , const std::function<double(double)>& fun
            , const std::function<double(double)>& funDeriv
            , int maxIter_ = 100
            , double toleranceThreshold_ = 1e-12
            , double zeroThreshold_ = 1e-20) {

            zeroThreshold = zeroThreshold_; 
            toleranceThreshold = toleranceThreshold_; 
            maxIter = maxIter_;
            x = x0;
            double xStep = 0;
            iter = 0;
            for (int i = 1; i <= maxIter; ++i) {

                iter += 1;
                fx = fun(x);
                dfx = funDeriv(x);
                if (std::abs(fx) < toleranceThreshold){
                    break;
                }
                
                if (std::abs(dfx) < zeroThreshold) {

                    x=std::numeric_limits<double>::quiet_NaN();
                    fx=std::numeric_limits<double>::quiet_NaN();
                    dfx=std::numeric_limits<double>::quiet_NaN();
                    break;
                }
                
                double xNew = x - fx / dfx;
                xStep = std::abs(xNew - x);
                x = xNew;

                if (xStep < toleranceThreshold) {break;}
                if (i==maxIter) break;
                
            }

        }

        NewtonRaphson(
            double x0
            , const std::function<std::pair<double,double>(double)>& funWithDeriv
            , int maxIter_ = 100
            , double toleranceThreshold_ = 1e-12
            , double zeroThreshold_ = 1e-20) {

            zeroThreshold = zeroThreshold_; 
            toleranceThreshold = toleranceThreshold_; 
            maxIter = maxIter_;
            x = x0;
            double xStep = 0;
            iter = 0;
            for (int i = 1; i <= maxIter; ++i) {

                iter += 1;
                std::pair<double,double> fwDerivAtX = funWithDeriv(x);
                fx = fwDerivAtX.first;
                dfx = fwDerivAtX.second;
                if (std::abs(fx) < toleranceThreshold){
                    break;
                }
                
                if (std::abs(dfx) < zeroThreshold) {

                    x=std::numeric_limits<double>::quiet_NaN();
                    fx=std::numeric_limits<double>::quiet_NaN();
                    dfx=std::numeric_limits<double>::quiet_NaN();
                    break;
                }
                
                double xNew = x - fx / dfx;
                xStep = std::abs(xNew - x);
                x = xNew;

                if (xStep < toleranceThreshold) {break;}
                if (i==maxIter) break;
                
            }

        }



    };
}