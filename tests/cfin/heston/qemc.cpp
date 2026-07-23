#include <numcpp/cfin/heston/qemc.hpp>
#include <numcpp/cfin/heston/cf.hpp>
#include <iostream>
#include <cassert>

double priceEuropeanOptionFromMatrix(const Eigen::MatrixXd& matrix, bool isCall, double strike) {

    size_t n = matrix.rows(), m = (matrix.cols()-1); 
    double totalPayoff = 0.0; 
    for (size_t i = 0; i<n;i++) {

        totalPayoff += isCall ? std::max(matrix(i,m)-strike,0.0) : std::max(strike-matrix(i,m),0.0);
    }
    return totalPayoff/double(n);
}

int main() {

    auto driftFun = [](double t) {return 0.05;};
    double T = 2.0;
    double S = 100.0;
    double F = S*std::exp(driftFun(2.0)*T);
    numcpp::cfin::Heston heston{0.05,0.05,0.45,.3,2.0};
    numcpp::cfin::HestonQuadraticExponentialSimulation sim(100.0,T,heston,driftFun, 2000, 500);
    double K = 100.0; 
    double x = std::log(F/K); 
    numcpp::gquad::GaussLaguerreQuadrature gl(128);

    std::cout << "Call MC price: " << priceEuropeanOptionFromMatrix(sim.xMatrix, true, K) << std::endl;
    std::cout << "Call semi-closed form price: " << heston.europeanNormalizedPriceGaussLaguerre(x,T,true,gl)*sqrt(F*K) << std::endl;
    std::cout << "Put MC price: " << priceEuropeanOptionFromMatrix(sim.xMatrix, false, K) << std::endl;
    std::cout << "Put semi-closed form price: " << heston.europeanNormalizedPriceGaussLaguerre(x,T,false,gl)*sqrt(F*K) << std::endl;

    return 0;
}