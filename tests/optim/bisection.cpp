#include <numcpp/optim/bisection.hpp>
#include <numcpp/optim/newtonraphson.hpp>
#include <cassert>


double f(double x) { return x * x - 2.0; };
double df(double x) { return 2.0 * x; };

void testBisectionVersusNewton() {

    numcpp::optim::NewtonRaphson nr(1.0,f,df); 
    numcpp::optim::Bisection bi(0.0, 2.0, f); 

    assert(std::abs(bi.x-nr.x) < 1e-4);

}

int main() {

    testBisectionVersusNewton();
    return 0; 
}