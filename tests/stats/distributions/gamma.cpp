#include <numcpp/stats/distributions/gamma.hpp>
#include <cassert>

inline bool isClose(double a, double b, double eps) { return std::abs(a-b)<eps;}

int main() {

    numcpp::stats::GammaDistribution g{2.0,2.0};
    
    assert(isClose(g.pdf(2.0),0.1839397206, 1e-9));
    assert(isClose(g.pdf(7.0),0.05284542099, 1e-9));
    assert(isClose(g.pdf(14.0),0.003191586879, 1e-9));

    assert(isClose(g.cdf(2.0),0.2642411177, 1e-9));
    assert(isClose(g.cdf(7.0),0.8641117746, 1e-9));
    assert(isClose(g.cdf(14.0),0.9927049443, 1e-9));
    return 0;
}