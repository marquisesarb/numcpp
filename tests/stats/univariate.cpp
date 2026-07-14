#include <numcpp/stats/univariate.hpp>
#include <cassert>

bool isClose(double a, double b, double tol = 1e-3) {return std::abs(a - b) <= tol;}

numcpp::stats::UnivariateStatistics univariateStatisticsObject() {

    Eigen::VectorXd data(50);
    data << 16, 14, 64, 88, 53, 3, 54, 96, 23, 88, 52,
        19, 92, 36, 94, 64, 69, 19, 92, 23, 96, 47,
        63, 17, 80, 59, 99, 7, 75, 23, 62, 96, 11, 78,
        70, 19, 32, 51, 44, 96, 90, 58, 64, 75, 20, 53,
        78, 95, 4, 82;

    return numcpp::stats::UnivariateStatistics{data};

}

numcpp::stats::UnivariateStatistics rollingMean11Result() {

    Eigen::VectorXd vec {{
        50.0909, 50.3636, 57.4545, 54.9091, 55.4545,
        56.4545, 62.4545, 59.2727, 58.9091, 58.9091,
        59.6364, 59.1818, 63.1818, 56.3636, 60.3636,
        57.1818, 60.3636, 54.7273, 59.8182, 53.5455,
        57.0909, 57.0909, 53.8182, 55.1818, 60.0000,
        54.4545, 52.0000, 47.6364, 51.0000, 52.9091,
        59.0000, 58.6364, 55.7273, 61.5455, 56.2727,
        54.7273, 60.0909, 65.8182, 61.5455
    }};

    return numcpp::stats::UnivariateStatistics{vec};

}

void testBasicStats() {

    numcpp::stats::UnivariateStatistics univobj = univariateStatisticsObject();

    assert(univobj.mean()==56.06);
    assert(isClose(univobj.variance(),925.8127,1e-4));
    assert(isClose(univobj.populationVariance(),907.2964,1e-4));
    assert(isClose(univobj.standardDeviation(),30.4272,1e-4));
    assert(isClose(univobj.populationStandardDeviation(),30.1214,1e-4));
    assert(isClose(univobj.skewness(),-0.2365,1e-4));
    assert(isClose(univobj.populationSkewness(),-0.2294,1e-4));
    assert(isClose(univobj.excessKurtosis(),-1.28317,1e-4));
    assert(isClose(univobj.populationExcessKurtosis(),-1.2760,1e-4));


    numcpp::stats::UnivariateStatistics rollMean11 = univobj.rollingMean(11);
    numcpp::stats::UnivariateStatistics rollMean11Result = rollingMean11Result();

    for (size_t i = 0; i<rollMean11.vector.size()-1; i++) {

        assert(isClose(rollMean11.vector(i), rollMean11Result.vector(i), 1e-4));
    }


    numcpp::reg::OLS ar = univobj.ar(1, true);

    assert(isClose(ar.betas(1), -0.3672, 1e-4)); 

    ar = univobj.ar(4, false);

    assert(isClose(ar.betas(0), -0.0745922, 1e-5)); 
    assert(isClose(ar.betas(1), 0.346834, 1e-5)); 
    assert(isClose(ar.betas(2), 0.34208, 1e-5)); 
    assert(isClose(ar.betas(3), 0.339985, 1e-5)); 

}



int main() {

    testBasicStats(); 
    return 0;
}

