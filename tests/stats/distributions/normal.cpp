#include <numcpp/stats/distributions/normal.hpp>
#include <cassert>

int main() {

    numcpp::stats::NormalDistribution n{}; 

    assert(std::abs(n.cdf(0.0) - 0.5) < 1e-6); 
    assert(std::abs(n.cdf(1.0) - 0.841344746) < 1e-6);
    assert(std::abs(n.cdf(-1.0) - 0.158655254) < 1e-6);
    assert(std::abs(n.cdf(2.0) - 0.977249868) < 1e-6);
    assert(std::abs(n.cdf(-2.0) - 0.022750132) < 1e-6);
    assert(std::abs(n.cdf(3.0) - 0.998650102) < 1e-6);

    assert(std::abs(n.invCdf(0.5) - 0.0) < 1e-6); 
    assert(std::abs(n.invCdf(.841344746) - 1.0) < 1e-6);
    assert(std::abs(n.invCdf(.158655254) + 1.0) < 1e-6);
    assert(std::abs(n.invCdf(0.977249868) - 2.0) < 1e-6);
    assert(std::abs(n.invCdf(0.022750132) + 2.0) < 1e-6);
    assert(std::abs(n.invCdf(0.998650102) - 3.0) < 1e-6);

    assert(std::abs(n.cf(0.0).real() - 1.0) < 1e-12);
    assert(std::abs(n.cf(0.0).imag()) < 1e-12);

    assert(std::abs(n.cf(0.1).real() - .995012479) < 1e-9);
    assert(std::abs(n.cf(0.1).imag()) < 1e-12);

    assert(std::abs(n.cf(1.0).real() - .6065306597) < 1e-9);
    assert(std::abs(n.cf(1.0).imag()) < 1e-12);

    assert(std::abs(n.pdf(0.0) - numcpp::stats::NormalDistribution::ONE_OVER_SQRT_TWO_PI) < 1e-6); // pdf(0) for standard normal

    n = numcpp::stats::NormalDistribution{2.0,3.0};

    assert(std::abs(n.cdf(2.0) - 0.5) < 1e-6);
    assert(std::abs(n.cdf(5.0) - 0.841344746) < 1e-6);
    assert(std::abs(n.cdf(-1.0) - 0.158655254) < 1e-6);

    assert(std::abs(n.invCdf(.5) - 2.0) < 1e-6);
    assert(std::abs(n.invCdf(.841344746) - 5.0) < 1e-6);
    assert(std::abs(n.invCdf(.158655254) + 1.0) < 1e-6);
   
    assert(std::abs(n.cf(0.0).real() - 1.0) < 1e-12);
    assert(std::abs(n.cf(0.0).imag()) < 1e-12);

    assert(std::abs(n.cf(0.1).real() - std::exp(-0.5 * 9.0 * 0.1 * 0.1) * std::cos(2.0 * 0.1)) < 1e-9);
    assert(std::abs(n.cf(0.1).imag() - std::exp(-0.5 * 9.0 * 0.1 * 0.1) * std::sin(2.0 * 0.1)) < 1e-9);

    assert(std::abs(n.pdf(2.0) - 1.0 / (3.0 * std::sqrt(2 * M_PI))) < 1e-6);


    return 0; 
}