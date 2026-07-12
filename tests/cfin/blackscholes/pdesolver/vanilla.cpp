#include <iostream>
#include <cassert>
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

    int N = 500; 
    int M = 100; 

    double absError = 1e9; 

    int MM = M; 

    double bsPrice = numcpp::cfin::blackScholesEuropeanNormalizedPrice(std::log(F/K), T, sigma,false)*std::sqrt(F*K)*std::exp(-r*T);

    while (MM < 600) {

        numcpp::cfin::BlackScholesPDESolverVanilla pdeEngine(S,K,T,mu,r,sigma,false,false,N,MM, useCrankNicolson,10);
        double intAbsError = std::abs(pdeEngine.price()-bsPrice);
        std::cout << "Absolute Error (M = " << MM << "): " << intAbsError << std::endl;
        absError = std::min(intAbsError,absError); 
        //assert(absError == intAbsError);
        MM += 100;
    }

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
    int M = 250; 
    
    std::cout << "Reference Price: " << ssvi.normalizedPrice(x, T, false)*std::sqrt(F*K)*std::exp(-r*T) << std::endl;
    std::cout << "Forward Price: " << F << std::endl;
    std::cout << "Spot Price: " << S*std::exp(-mu*T) << std::endl;

    numcpp::cfin::BlackScholesPDESolverVanilla pdeEngine(S,K,T,true,false,driftFun,localVolFun,discountRateFun,N,M, useCrankNicolson,5);

  
}


int main() {

    
    europeanFlatSigma(false);
    europeanFlatSigma(true);
    //europeanLocalVolatility(false);



}