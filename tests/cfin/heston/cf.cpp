#include <numcpp/cfin/heston/cf.hpp>
#include <cassert>

bool isClose(double a, double b, double tol = 1e-3) {return std::abs(a - b) <= tol;}

int main() {

    double S = 100.0;    
    double K = 100.0;    
    double r = 0.05;   
    double q = 0.01;    
    double T = 1.5;     
    double F = S * std::exp((r - q) * T);

    numcpp::gquad::GaussLaguerreQuadrature gl(64);
    numcpp::cfin::Heston heston{0.05,0.05,0.45,.3,2.0};
    double expectedHestonPrice = 13.2561; 
    double x = std::log(F/K);
    double modelHestonPrice = heston.europeanNormalizedPriceGaussLaguerre(x,T,true,gl)*std::exp(-r*T)*sqrt(F*K);
    assert(isClose(modelHestonPrice, expectedHestonPrice, 1e-4));
    return 0;
}