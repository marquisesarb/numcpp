#pragma once 
#include <cmath>
#include <numcpp/stats/distributions/interface.hpp>

namespace numcpp::stats {

    struct UniformDistribution final : public ContinuousUnivariateProbabilityDistribution {

        double a=0.0, b=1.0; 

        UniformDistribution(double a_, double b_): a(a_), b(b_) {}
        UniformDistribution(): a(0.0), b(1.0) {}

        double pdf(double x) const override {return (x>=a and x<=b) ? 1.0/(b-a) : 0.0;}
        double cdf(double x) const override {return (x<a) ? 0.0 : (x>b) ? 1.0 : (x-a)/(b-a);}  

        std::complex<double> cf(double t) const override {
            using cd = std::complex<double>;
            const cd i(0.0, 1.0);
            if (t == 0.0) { return cd(1.0, 0.0);}
            return (std::exp(i * t * b) - std::exp(i * t * a)) / (i * t * (b - a));
        }

        double sample(std::mt19937& gen) const override {
            std::uniform_real_distribution<> dist(a, b);
            return dist(gen);
        }

        inline double invCdf(double p) const {return a + p*(b-a);}

        double mean() const {return .5*(a+b);}
        double variance() const {return (b-a)*(b-a)/12.0;}
        double skewness() const {return 0.0;}
        double kurtosis() const {return -6.0/5.0;}

    };


}