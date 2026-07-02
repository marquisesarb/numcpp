#include <numcpp/stats/distributions/chisq.hpp>
#include <cassert>

inline bool isClose(double a, double b, double eps) { return std::abs(a-b)<eps;}

int main() {

    numcpp::stats::ChiSquaredDistribution g{2.0};
    
    assert(isClose(g.pdf(2.0),0.1839397206, 1e-9));
    assert(isClose(g.pdf(7.0),0.01509869171, 1e-9));
    assert(isClose(g.pdf(14.0),0.0004559409828, 1e-9));

    assert(isClose(g.pdf(2.0),0.6321205588, 1e-9));
    assert(isClose(g.pdf(7.0),0.9698026166, 1e-9));
    assert(isClose(g.pdf(14.0),0.999088118, 1e-9));
    return 0;
}