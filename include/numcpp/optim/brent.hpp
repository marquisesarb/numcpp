#pragma once 
#include <cmath>
#include <functional>

namespace numcpp::optim {

    struct Brent {

        double x;
        double fx;
        int iter;
        int maxIter; 
        double toleranceThreshold; 
        bool converged;

        Brent(
            double ax, 
            double bx, 
            const std::function<double(double)>& fun, 
            int maxIter_ = 100, 
            double toleranceThreshold_ = 1e-12) {

            maxIter = maxIter_; 
            toleranceThreshold = toleranceThreshold_;
            iter = 0; 

            const double phi = 0.5 * (3.0 - std::sqrt(5.0));

            double a = ax;
            double b = bx;

            x = a + phi * (b - a);
            fx = fun(x);

            double w = x;
            double v = x;
            double fw = fx;
            double fv = fx;

            double d = 0.0;  
            double e = 0.0;  

            while (true) {

                double m = 0.5 * (a + b);
                double tol1 = toleranceThreshold * std::abs(x) + 1e-12;
                double tol2 = 2.0 * tol1;

                if (std::abs(x - m) <= tol2 - 0.5 * (b - a))  {converged = true; break;}

                bool parabolicStepAccepted = false;

                if (std::abs(e) > tol1) {
                    
                    double r = (x - w) * (fx - fv);
                    double q = (x - v) * (fx - fw);
                    double p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);

                    if (q > 0.0) p = -p;
                    q = std::abs(q);

                    if (std::abs(p) < std::abs(0.5 * q * e) &&
                        p > q * (a - x) &&
                        p < q * (b - x)) {

                        d = p / q;
                        parabolicStepAccepted = true;
                    }
                }

                if (!parabolicStepAccepted) {
                    e = (x < m) ? (b - x) : (a - x);
                    d = phi * e;
                }

                double u = x + ((std::abs(d) >= tol1) ? d : (d > 0 ? tol1 : -tol1));
                double fu = fun(u);

                if (fu <= fx) {
                    if (u < x) b = x;
                    else       a = x;

                    v = w; fv = fw;
                    w = x; fw = fx;
                    x = u; fx = fu;
                } else {
                    if (u < x) a = u;
                    else       b = u;

                    if (fu <= fw || w == x) {
                        v = w; fv = fw;
                        w = u; fw = fu;
                    } else if (fu <= fv || v == x || v == w) {
                        v = u; fv = fu;
                    }
                }

                if (iter == maxIter) {converged = false; break;}

                e = d;
                iter++;
            }            
        }
    };

}