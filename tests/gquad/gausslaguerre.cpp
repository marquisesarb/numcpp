#include <numcpp/gquad/gausslaguerre.hpp>
#include <cassert>

void testGaussLaguerre() {

    int points = 20;
    numcpp::gquad::GaussLaguerreQuadrature glq(points);

    assert(glq.roots.size() == points);
    assert(glq.weights.size() == points);

    auto f = [](double x) { return x; }; 
    double result = glq.integrate(f);

    assert(std::abs(result - 1.0) < 1e-6);
}

void testGaussLaguerreHighDimension() {

    int points = 300;
     numcpp::gquad::GaussLaguerreQuadrature glq(points);

    assert(glq.roots.size() == points);
    assert(glq.weights.size() == points);

    auto f = [](double x) { return x; }; 
    double result = glq.integrate(f);

    assert(std::abs(result - 1.0) < 1e-6);
}

int main() {
    testGaussLaguerre();
    testGaussLaguerreHighDimension();
    return 0;
}