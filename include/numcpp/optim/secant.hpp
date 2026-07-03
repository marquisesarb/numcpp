#pragma once 
#include <cmath>
#include <functional>

namespace numcpp::optim {

    struct Secant {

        double x;
        double fx;
        int iter;
        int maxIter; 
        double toleranceThreshold; 
        bool converged;

        Secant(
            double x0,
            double x1,
            const std::function<double(double)>& f,
            int maxIter_ = 100,
            double toleranceThreshold_ = 1e-12) {

            maxIter = maxIter_; 
            toleranceThreshold = toleranceThreshold_;
            iter = 0; 

            fx = f(x1);
            double fx_1 = f(x0);
            x = x1;
            double x_1 = x0;
            
            while (true) {
                

                double ix_1 = x_1; 
                x_1 = x;
                x = (ix_1*fx-x*fx_1)/(fx-fx_1);
                fx_1 = fx;
                fx = f(x);
                iter++;
            
                if (std::abs(fx) < toleranceThreshold) {converged = true; break;}
                if (std::abs(x1 - x0) < toleranceThreshold) {converged = true; break;}
                if (iter==maxIter) {{converged = false; break;}}
            }
        }
    };
}