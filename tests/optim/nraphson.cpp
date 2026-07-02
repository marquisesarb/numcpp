#include <numcpp/optim/newtonraphson.hpp>
#include <cassert>

double testFunctionNR(double x) {
    return x * x - 4;  
}

double testFunctionDerivativeNR(double x) {
    return 2 * x;
}

std::pair<double,double> testFunctionWithDerivativeNR(double x) {

    return {x * x - 4, 2*x};
}

int main() {


    numcpp::optim::NewtonRaphson nr(1.0, testFunctionNR, testFunctionDerivativeNR);
    assert(std::abs(nr.x - 2.0) < 1e-4 || std::abs(nr.x+ 2.0) < 1e-4);

    numcpp::optim::NewtonRaphson nr2(1.0, testFunctionWithDerivativeNR);
    assert(std::abs(nr2.x - 2.0) < 1e-4 || std::abs(nr2.x+ 2.0) < 1e-4);

    return 0; 
}