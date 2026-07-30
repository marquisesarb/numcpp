#include "Eigen/Core"
#include <numcpp/stats/univariate.hpp>
#include <limits>

namespace numcpp::stats {

    double sum(const Eigen::VectorXd& data) {return data.sum();}

    Eigen::VectorXd centered(const Eigen::VectorXd& data) {

        return data.array() - mean(data);
    }

    Eigen::VectorXd normalized(const Eigen::VectorXd& data) {

        return (data.array() - mean(data))/standardDeviation(data);
    }

    double mean(const Eigen::VectorXd& data) { return data.size()>0 ? data.mean() : std::numeric_limits<double>::quiet_NaN();}

    double variance(const Eigen::VectorXd& data) {

        if (data.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        double M2 =0.0; 
        double mean = 0.0, oldmean =0.0; 
        for (size_t i =0; i<data.size();i++ ) {
            oldmean = mean; 
            mean += (data(i)-oldmean)/ double(i+1); 
            M2 += (data(i)-oldmean)*(data(i)-mean);
        }
        return M2/double(data.size()-1);
    } 

    double populationVariance(const Eigen::VectorXd& data) { 

        if (data.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        return variance(data)*double(data.size()-1)/double(data.size());
    }

    double standardDeviation(const Eigen::VectorXd& data) {

        if (data.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(variance(data));
    }

    double populationStandardDeviation(const Eigen::VectorXd& data) {
        
        if (data.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(populationVariance(data));
    }

    double populationSkewness(const Eigen::VectorXd& data) {

        if (data.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        Eigen::ArrayXd centered_ = data.array() - mean(data);
        double m2 = centered_.square().mean();
        return (m2==0.0) ? 0.0 : centered_.pow(3).mean() / std::pow(m2, 1.5); 
    }

    double skewness(const Eigen::VectorXd& data) {

        int n = data.size();
        if (n == 0 ||  n < 3) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(n * (n - 1.0)) / (n - 2.0)*populationSkewness(data);
    }

    double populationExcessKurtosis(const Eigen::VectorXd& data) {

        int n = data.size();
        if (n == 0) return std::numeric_limits<double>::quiet_NaN();
        Eigen::ArrayXd centered_ = data.array() - mean(data);
        double m2 = centered_.square().mean();
        if (m2 == 0.0) return -3.0;
        return (m2==0.0) ? -3.0 : centered_.pow(4).mean() / (m2 * m2) - 3.0;  
    }

    double excessKurtosis(const Eigen::VectorXd& data) {

        int n = data.size();
        if (n == 0 ||  n < 4) return std::numeric_limits<double>::quiet_NaN();
        double n1 = (n - 1.0);
        double n2 = (n - 2.0);
        double n3 = (n - 3.0);
        return ((n - 1.0) / ((n - 2.0) * (n - 3.0))) *
            ((n + 1.0) * populationExcessKurtosis(data) + 6.0);
    }

    double autoCorrelation(const Eigen::VectorXd& data, size_t lag) {

        if (lag<0 or lag > data.size()) return std::numeric_limits<double>::quiet_NaN();
        int n = data.size()-lag;
        Eigen::VectorXd x = data.segment(lag, n); 
        Eigen::VectorXd y = data.segment(0, n); 

        double sumxy =0.0; 
        double sumx =0.0;
        double sumy =0.0; 
        double sumysq =0.0; 
        double sumxsq =0.0; 


        for (size_t i=0; i<n; i++) {
            sumxy += x(i)*y(i);
            sumxsq += x(i)*x(i);
            sumx += x(i); 
            sumysq += y(i)*y(i);
            sumy += y(i); 
        }

        return (n*sumxy - sumx*sumy)/std::sqrt((n*sumxsq-sumx*sumx)*(n*sumysq-sumy*sumy));
        

    }

    numcpp::reg::OLS ar(const Eigen::VectorXd& data, size_t lag, bool partialCorrelation = false, bool intercept = false) {

        size_t xSize = data.size()-lag;
        Eigen::MatrixXd X = Eigen::MatrixXd::Zero(xSize, lag);
        Eigen::VectorXd seg = data.segment(lag,xSize);
        Eigen::VectorXd Y = partialCorrelation ? normalized(seg) : seg;

        for (size_t i =0; i<lag; i++) {

            seg = data.segment(i,xSize);
            X.col(lag-1-i) = partialCorrelation ? normalized(seg) : seg;
        }

        return numcpp::reg::OLS(Y, X, intercept);
    }

    Eigen::VectorXd rollingTool(const Eigen::VectorXd& data, size_t window, const std::function<double(const Eigen::VectorXd&)> lambda_) {

        Eigen::VectorXd newVector(data.size()-window+1);

        for (size_t i = 0; i<data.size()-window+1; i++) {

            newVector(i) = lambda_(data.segment(i,window));
        }

        return newVector;
    }

    Eigen::VectorXd rollingMean(const Eigen::VectorXd& data, size_t window) {return rollingTool(data, window, [](const Eigen::VectorXd& u) {return u.mean();});}



}