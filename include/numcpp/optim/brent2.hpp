#pragma once 
#include <functional>
#include <cmath>

namespace numcpp::optim {


    struct Brent {

        double x;
        double fx;
        int iter;
        int maxIter; 
        double toleranceThreshold; 
        bool converged;

        Brent(
            double a,
            double b,
            const std::function<double(double)>& f,
            int maxIter_ = 100,
            double toleranceThreshold_ = 1e-12) {

            maxIter = maxIter_; 
            toleranceThreshold = toleranceThreshold_;
            iter = 0; 
            double c = b; 
            double fa = f(a); 
            double fb = f(b); 
            double fc = f(c); 
            

            while (true) {

                ;
            }


        }
    };

}