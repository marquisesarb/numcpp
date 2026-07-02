#pragma once 
#include <complex>
#include <limits>
#include <cmath>
#include <numcpp/stats/distributions/chisq.hpp>
#include <numcpp/optim/newtonraphson.hpp>


namespace numcpp::stats {

    inline double logbeta(double a,double b) {return std::lgamma(a + b)- std::lgamma(a)- std::lgamma(b);}

    inline double incompleteBetaContinuedFraction(double a, double b, double x) {

        double qab = a + b;
        double qap = a + 1.0;
        double qam = a - 1.0;

        double c = 1.0;
        double d = 1.0 - qab * x / qap;
        if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
        d = 1.0 / d;
        double h = d;

        for (int m = 1; m <= 200; ++m)
        {
            int m2 = 2 * m;

            // even step
            double aa = m * (b - m) * x /
                        ((qam + m2) * (a + m2));

            d = 1.0 + aa * d;
            if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
            c = 1.0 + aa / c;
            if (std::fabs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
            d = 1.0 / d;
            h *= d * c;

            // odd step
            aa = -(a + m) * (qab + m) * x /
                ((a + m2) * (qap + m2));

            d = 1.0 + aa * d;
            if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
            c = 1.0 + aa / c;
            if (std::fabs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
            d = 1.0 / d;
            double del = d * c;
            h *= del;

            if (std::fabs(del - 1.0) < 1e-14)
                break;
        }

        return h;
    }

    inline double incompleteBeta(double a, double b, double x) {
        
        if (x == 0.0) return 0.0;
        if (x == 1.0) return 1.0;

        bool flip = false;
        if (x > (a + 1.0) / (a + b + 2.0))
        {
            flip = true;
            std::swap(a, b);
            x = 1.0 - x;
        }

        double front = std::exp(
            a * std::log(x) +
            b * std::log(1.0 - x) +
            logbeta(a,b)
        ) / a;

        double cf = incompleteBetaContinuedFraction(a, b, x);
        double result = front * cf;

        return flip ? 1.0 - result : result;
    }

    inline double cornishFisherStudentInverseCdfExpansion(double p, double nu) {

        double z = acklamStandardGaussianInverseCdf(p);

        double z2 = z*z;
        double z3 = z2*z;
        double z5 = z3*z2;

        return z + (z3 + z)/(4.0*nu)
                    + (5.0*z5 + 16.0*z3 + 3.0*z)/(96.0*nu*nu);
    }
    
    struct TStudentDistribution final : public ContinuousUnivariateProbabilityDistribution {

        static constexpr double PI = 3.14159265358979323846;
        double nu;

        TStudentDistribution(double nu_): nu(nu_) {}

        double pdf(double x) const override { return std::pow(1.0+x*x/nu,-.5*(nu+1.0))*std::tgamma(.5*(nu+1.0))/(std::sqrt(nu*PI)*std::tgamma(.5*nu));}
        
        double cdf(double x) const override { return x>=0.0 ? 1.0-.5*incompleteBeta(nu/2.0,.5,nu/(nu+x*x)) : .5*incompleteBeta(nu/2.0,.5,nu/(nu+x*x));}

        double sample(std::mt19937& gen) const override {return NormalDistribution{}.sample(gen)/std::sqrt(ChiSquaredDistribution{nu}.sample(gen)/nu);}

        std::complex<double> cf(double t) const override {
            // must have the modified bessel function, see : https://en.wikipedia.org/wiki/Student%27s_t-distribution
            return std::numeric_limits<std::complex<double>>::quiet_NaN();
        }

        double invCdf(double p) const {

            double z0 = cornishFisherStudentInverseCdfExpansion(p, nu);
            std::function<double(double)> target = [*this,p](double x) { return p-cdf(x); };
            std::function<double(double)> targetDeriv = [*this,p](double x) { return -pdf(x); };
            optim::NewtonRaphson result(z0, target, targetDeriv,20,1e-12);
            return result.x;
            //return std::numeric_limits<double>::quiet_NaN();
        }

    };

}