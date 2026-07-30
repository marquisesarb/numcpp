#pragma once 
#include <functional>
#include <Eigen/Dense>
#include <numcpp/regression/ols.hpp>

namespace numcpp::stats {

    double sum(const Eigen::VectorXd& data);
    double mean(const Eigen::VectorXd& data);
    double variance(const Eigen::VectorXd& data);
    double populationVariance(const Eigen::VectorXd& data);
    double standardDeviation(const Eigen::VectorXd& data);
    double populationStandardDeviation(const Eigen::VectorXd& data);
    double populationSkewness(const Eigen::VectorXd& data);
    double skewness(const Eigen::VectorXd& data);
    double populationExcessKurtosis(const Eigen::VectorXd& data);
    double excessKurtosis(const Eigen::VectorXd& data);
    double autoCorrelation(const Eigen::VectorXd& data, size_t lag);
    Eigen::VectorXd normalized(const Eigen::VectorXd& data);
    Eigen::VectorXd centered(const Eigen::VectorXd& data);
    Eigen::VectorXd rollingTool(size_t window, const std::function<double(const Eigen::VectorXd&)> lambda_); 
    Eigen::VectorXd rollingMean(const Eigen::VectorXd& data, size_t window);
    numcpp::reg::OLS ar(const Eigen::VectorXd& data, size_t lag,  bool partialCorrelation, bool intercept);

}