#pragma once 
#include <cmath>
#include <functional>
#include <limits>

namespace numcpp::optim {

    struct Bisection {

        double x;
        double fx;
        int iter;
        int maxIter; 
        double toleranceThreshold; 
        bool converged;

        Bisection(
            double a,
            double b,
            const std::function<double(double)>& fun,
            int maxIter_ = 100,
            double toleranceThreshold_ = 1e-12) {

            maxIter = maxIter_; 
            toleranceThreshold = toleranceThreshold_;
            iter = 0; 

            double fa = fun(a);
            double fb = fun(b);

            if (fa * fb > 0) {

                x = std::numeric_limits<double>::quiet_NaN(); 
                fx = std::numeric_limits<double>::quiet_NaN(); 
                converged = false;

            } else {

                x = a;
                fx = fa;

                while (true) {
                    

                    x = 0.5 * (a + b);
                    fx = fun(x);
                    iter++;

                    if (std::abs(fx) < toleranceThreshold) {converged = true; break;}
                    if (std::abs(b - a) < toleranceThreshold) {converged = true; break;}

                    if (fa * fx < 0) {
                        b = x;
                        fb = fx;
                    } else {
                        a = x;
                        fa = fx;
                    }

                    if (iter == maxIter){converged = false; break;}
                }
            } 
        }
    };
}