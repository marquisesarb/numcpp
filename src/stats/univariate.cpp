#include <numcpp/stats/univariate.hpp>
#include <limits>

namespace numcpp::stats {

    double UnivariateStatistics::sum() const {return vector.sum();}

    UnivariateStatistics UnivariateStatistics::centered() const {

        return UnivariateStatistics{vector.array() - mean()};
    }

    UnivariateStatistics UnivariateStatistics::normalized() const {

        Eigen::VectorXd array = vector.array() - mean(); 
        Eigen::VectorXd normArray = array/standardDeviation();
        return UnivariateStatistics{array};
    }

    UnivariateStatistics UnivariateStatistics::segment(size_t start, size_t n) const {

        if (start<0 or start > vector.size()-1) return UnivariateStatistics{};
        return UnivariateStatistics{vector.segment(start, n)};
    }

    double UnivariateStatistics::mean() const { return vector.size()>0 ? vector.mean() : std::numeric_limits<double>::quiet_NaN();}

    double UnivariateStatistics::variance() const {

        if (vector.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        double M2 =0.0; 
        double mean = 0.0, oldmean =0.0; 
        for (size_t i =0; i<vector.size();i++ ) {
            oldmean = mean; 
            mean += (vector(i)-oldmean)/ double(i+1); 
            M2 += (vector(i)-oldmean)*(vector(i)-mean);
        }
        return M2/double(vector.size()-1);
    } 

    double UnivariateStatistics::populationVariance() const { 

        if (vector.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        return variance()*double(vector.size()-1)/double(vector.size());
    }

    double UnivariateStatistics::standardDeviation() const {

        if (vector.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(variance());
    }

    double UnivariateStatistics::populationStandardDeviation() const {
        
        if (vector.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(populationVariance());
    }

    double UnivariateStatistics::populationSkewness() const {

        if (vector.size()<=0) return std::numeric_limits<double>::quiet_NaN();
        Eigen::ArrayXd centered_ = vector.array() - mean();
        double m2 = centered_.square().mean();
        return (m2==0.0) ? 0.0 : centered_.pow(3).mean() / std::pow(m2, 1.5); 
    }

    double UnivariateStatistics::skewness() const {

        int n = vector.size();
        if (n == 0 ||  n < 3) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(n * (n - 1.0)) / (n - 2.0)*populationSkewness();
    }

    double UnivariateStatistics::populationExcessKurtosis() const {

        int n = vector.size();
        if (n == 0) return std::numeric_limits<double>::quiet_NaN();
        Eigen::ArrayXd centered_ = vector.array() - mean();
        double m2 = centered_.square().mean();
        if (m2 == 0.0) return -3.0;
        return (m2==0.0) ? -3.0 : centered_.pow(4).mean() / (m2 * m2) - 3.0;  
    }

    double UnivariateStatistics::excessKurtosis() const {

        int n = vector.size();
        if (n == 0 ||  n < 4) return std::numeric_limits<double>::quiet_NaN();
        double n1 = (n - 1.0);
        double n2 = (n - 2.0);
        double n3 = (n - 3.0);
        return ((n - 1.0) / ((n - 2.0) * (n - 3.0))) *
            ((n + 1.0) * populationExcessKurtosis() + 6.0);
    }

    numcpp::reg::OLS UnivariateStatistics::ar(size_t lag, bool intercept = false) const {

        size_t xSize = vector.size()-lag;
        Eigen::MatrixXd X = Eigen::MatrixXd::Zero(xSize, lag);
        Eigen::VectorXd Y = vector.segment(lag,xSize);

        for (size_t i =0; i<lag; i++) {

            X.col(lag-1-i) = vector.segment(i,xSize);
        }

        return numcpp::reg::OLS(Y, X, intercept);
    }


    UnivariateStatistics UnivariateStatistics::rollingTool(size_t window, const std::function<double(const UnivariateStatistics&)> lambda_) const {

        Eigen::VectorXd newVector(vector.size()-window+1);

        for (size_t i = 0; i<vector.size()-window+1; i++) {

            newVector(i) = lambda_(segment(i,window));
        }

        

        return UnivariateStatistics{newVector};

    }

    UnivariateStatistics UnivariateStatistics::rollingMean(size_t window) const {return rollingTool(window, [*this](const UnivariateStatistics& u) {return u.mean();});}



}