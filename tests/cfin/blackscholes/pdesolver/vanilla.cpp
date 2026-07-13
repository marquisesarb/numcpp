#include <iostream>
#include <cassert>
#include <chrono>
#include <numcpp/cfin/blackscholes/pdesolver/vanilla.hpp>
#include <numcpp/cfin/blackscholes/cf.hpp>
#include <numcpp/cfin/blackscholes/baw.hpp>
#include <numcpp/cfin/ssvi.hpp>


void europeanFlatSigma(bool useCrankNicolson) {

    double S = 100.0; 
    double r = 0.001;
    double mu = r - .11;
    double T = 1.0;
    double F = S*std::exp(mu*T); 
    double K = 100.0;
    double sigma = 0.16; 

    int N = 200; 
    int M = 125; 

    double bsPrice = numcpp::cfin::blackScholesEuropeanNormalizedPrice(std::log(F/K), T, sigma,false)*std::sqrt(F*K)*std::exp(-r*T);

    std::cout << "Closed Form Price: " << bsPrice << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla pdeEngine(S,K,T,mu,r,sigma,false,false,N,M, useCrankNicolson,10);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "BS PDE Price: " << pdeEngine.price() << std::endl;
    std::cout << "Time taken for solving BS PDE: " << ms_int.count()/1000.0 << "s\n";

}

void europeanLocalVolatility(bool useCrankNicolson) {

    numcpp::cfin::SSVIPowerLawFlatATM ssvi(-0.5,0.46,0.77,1.0); 

    double mu = .05;
    double r = 0.01; 

    auto driftFun = [mu](double t) {return mu;};
    auto discountRateFun = [r](double t) {return r;};
    auto localVolFun = [ssvi](double k, double t) {return std::sqrt(ssvi.localVariance(k, t));};

    double T = 2.0;
    double S = 100.0;
    double F = S*std::exp(mu*T);
    double K = 100.0; 
    double x = std::log(F/K); 
    int N = 400; 
    int M = 200; 
    
    std::cout << "Closed Form Price: " << ssvi.normalizedPrice(x, T, false)*std::sqrt(F*K)*std::exp(-r*T) << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla pdeEngine(S,K,T,false,false,driftFun,localVolFun,discountRateFun,N,M, useCrankNicolson,5);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local Vol PDE Price: " << pdeEngine.price() << std::endl;
    std::cout << "Time taken for solving Local-Vol PDE: " << ms_int.count() << "ms\n";
  
}


int main() {

    europeanFlatSigma(false);

    europeanLocalVolatility(false);
    //europeanLocalVolatility(false);



}