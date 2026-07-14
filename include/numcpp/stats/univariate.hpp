#pragma once 
#include <functional>
#include <Eigen/Dense>
#include <numcpp/regression/ols.hpp>


namespace numcpp::stats {

    struct UnivariateStatistics {

        Eigen::VectorXd vector;

        double sum() const;
        double mean() const;
        double variance() const;
        double populationVariance() const;
        double standardDeviation() const;
        double populationStandardDeviation() const;
        double populationSkewness() const;
        double skewness() const;
        double populationExcessKurtosis() const;
        double excessKurtosis() const;
        double autoCorrelation(size_t lag) const;
        UnivariateStatistics normalized() const;
        UnivariateStatistics centered() const;
        UnivariateStatistics segment(size_t start, size_t n) const;
        UnivariateStatistics rollingTool(size_t window, const std::function<double(const UnivariateStatistics&)> lambda_) const; 
        UnivariateStatistics rollingMean(size_t window) const;
        numcpp::reg::OLS ar(size_t lag,  bool intercept) const;
         
    };

}