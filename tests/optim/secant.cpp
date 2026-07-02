#include <numcpp/optim/secant.hpp>
#include <numcpp/optim/newtonraphson.hpp>
#include <cassert>

double f(double x) { return x * x - 612.0; };
double df(double x) { return 2.0*x; };

int main() {

    numcpp::optim::Secant sec(10,30,f,5,1e-12); 
    numcpp::optim::NewtonRaphson nr(10,f,df,10,1e-12,1e-20);
    assert(std::abs(nr.x-sec.x)<1e-4);
    assert(std::abs(nr.fx-sec.fx)<1e-4);
}