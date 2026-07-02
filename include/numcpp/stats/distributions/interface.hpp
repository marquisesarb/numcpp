#pragma once 
#include <random>
#include <complex>
#include <Eigen/Dense>

namespace numcpp::stats {

    struct ProbabilityDistribution {

        virtual ~ProbabilityDistribution() = default; 

        virtual double cdf(double x) const = 0;  
        virtual std::complex<double> cf(double t) const = 0; 

    };

    struct UnivariateProbabilityDistribution : public ProbabilityDistribution {

        virtual ~UnivariateProbabilityDistribution() = default; 

        virtual double cdf(double x) const override = 0;  
        virtual std::complex<double> cf(double t) const override = 0; 
        virtual double sample(std::mt19937& gen) const = 0;  

        Eigen::VectorXd samples(std::mt19937& gen, int n) const {
            return Eigen::VectorXd::NullaryExpr(n, [&]() {return sample(gen);});
        }

        Eigen::MatrixXd samples(std::mt19937& gen, int n, int m) const {
            return Eigen::MatrixXd::NullaryExpr(n, m, [&]() {return sample(gen);});
        }
    };


    struct ContinuousUnivariateProbabilityDistribution : public UnivariateProbabilityDistribution {

        virtual ~ContinuousUnivariateProbabilityDistribution() = default; 

        virtual double pdf(double x) const = 0;
        virtual double cdf(double x) const override = 0;  
        virtual std::complex<double> cf(double t) const override = 0; 
        virtual double sample(std::mt19937& gen) const override = 0; 

    };

    struct DiscreteUnivariateProbabilityDistribution : public UnivariateProbabilityDistribution {

        virtual ~DiscreteUnivariateProbabilityDistribution() = default; 

        virtual double pmf(double x) const = 0;
        virtual double cdf(double x) const override = 0;  
        virtual std::complex<double> cf(double t) const override = 0; 
        virtual double sample(std::mt19937& gen) const override = 0; 

    };


}