#include <numcpp/stats/distributions/tstudent.hpp>
#include <cassert>

inline bool isClose(double a, double b, double eps) { return std::abs(a-b)<eps;}

int main() {

    numcpp::stats::TStudentDistribution g{2.0};
    
    assert(isClose(g.pdf(-8.0),0.001865022591, 1e-9));
    assert(isClose(g.pdf(-2.0),0.06804138174, 1e-9));
    assert(isClose(g.pdf(.0),0.3535533906, 1e-9));
    assert(isClose(g.pdf(3.0),0.02741012223, 1e-9));
    assert(isClose(g.pdf(9.0),0.001322460964, 1e-9));

    assert(isClose(g.cdf(-8.0),0.007634036083, 1e-9));
    assert(isClose(g.cdf(-2.0),0.09175170954, 1e-9));
    assert(isClose(g.cdf(.0),0.5, 1e-9));
    assert(isClose(g.cdf(3.0),0.9522670169, 1e-9));
    assert(isClose(g.cdf(9.0),0.99393917, 1e-9));

    assert(isClose(g.invCdf(0.0000001),-2236.067642, 1e-3));
    assert(isClose(g.invCdf(.01),-6.964556734, 1e-9));
    assert(isClose(g.invCdf(.2),-1.060660172, 1e-9));
    assert(isClose(g.invCdf(0.9999),70.70007107, 1e-4));
    assert(isClose(g.invCdf(0.9999999),2236.067642, 1e-3));

    return 0;
}
