#include <numcpp/cfin/blackscholes/mc.hpp>
#include <numcpp/cfin/blackscholes/cf.hpp>
#include <iostream>

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
    double T = 2.0, sigma = 0.2;
    double S = 100.0;
    double mu = .05;
    double F = S*std::exp(mu*T);
    numcpp::cfin::BlackScholesEulerSimulation sim(100.0,T,sigma,mu,1000,500);
    double K = 100.0; 
    double x = std::log(F/K); 
    std::cout << "Call closed form price: " << numcpp::cfin::blackScholesEuropeanNormalizedPrice(x, T, sigma, true)*std::sqrt(K*F) << std::endl;
    std::cout << "Call monte carlo price: " << priceEuropeanOptionFromMatrix(sim.xMatrix, true, K) << std::endl;
    std::cout << "Put closed form price: " << numcpp::cfin::blackScholesEuropeanNormalizedPrice(x, T, sigma, false)*std::sqrt(K*F) << std::endl;
    std::cout << "Put monte carlo price: " << priceEuropeanOptionFromMatrix(sim.xMatrix, false, K) << std::endl;
    return 0;
}