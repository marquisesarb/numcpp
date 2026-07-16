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
    int M = 75; 

    double bsPrice = numcpp::cfin::blackScholesEuropeanNormalizedPrice(std::log(F/K), T, sigma,false)*std::sqrt(F*K)*std::exp(-r*T);

    
    auto t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla pdeEngine(S,K,T,mu,r,sigma,false,false,N,M, useCrankNicolson,10);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    std::cout << "Closed Form Price: " << bsPrice << std::endl;
    std::cout << "BS PDE Price: " << pdeEngine.price() << std::endl;
    std::cout << "Time taken for solving BS PDE: " << ms_int.count() << "ms\n";

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
    int N = 200; 
    int M = 125; 
    
    std::cout << "Closed Form Price: " << ssvi.normalizedPrice(x, T, false)*std::sqrt(F*K)*std::exp(-r*T) << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla pdeEngine(S,K,T,false,false,driftFun,localVolFun,discountRateFun,N,M, useCrankNicolson,5);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local Vol PDE Price: " << pdeEngine.price() << std::endl;
    std::cout << "Time taken for solving Local-Vol PDE: " << ms_int.count() << "ms\n";
  
}

void functionMdFile() {

     numcpp::cfin::SSVIPowerLawFlatATM ssvi(-0.5,0.46,0.77,1.0); 

    double mu = -.11;
    double r = 0.07; 

    auto driftFun = [mu](double t) {return mu;};
    auto discountRateFun = [r](double t) {return r;};
    auto localVolFun = [ssvi](double k, double t) {return std::sqrt(ssvi.localVariance(k, t));};

    double T = 2.0;
    double S = 100.0;
    double F = S*std::exp(mu*T);
    double K = 100.0; 
    double x = std::log(F/K); 
    int N = 200; 
    int M = 125; 
    double sigma = std::sqrt(ssvi.totalVariance(x,T)/T);

    std::cout << "Closed Form European Put Price: " << ssvi.normalizedPrice(x, T, false)*std::sqrt(F*K)*std::exp(-r*T) << std::endl;
    // Closed Form Price: 48.358 

    auto t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanLocalVol(S,K,T,false,false,driftFun,localVolFun,discountRateFun,N,M, false,10);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local Vol PDE European Put Price: " << europeanLocalVol.price() << " (" << msEuropeanLocalVol.count() << " milliseconds)" << std::endl;
    //

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanFlatVol(S,K,T,mu,r,sigma,false,false,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE European Put Price: " << europeanFlatVol.price() << " (" << msEuropeanFlatVol.count() << " milliseconds)" << std::endl;
    //

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanLocalVol(S,K,T,false,true,driftFun,localVolFun,discountRateFun,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local-Vol PDE American Put Price: " << americanLocalVol.price() << " (" << msAmericanLocalVol.count() << " milliseconds)" << std::endl;
    //

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanFlatVol(S,K,T,mu,r,sigma,false,true,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE American Put Price: " << americanFlatVol.price() << " (" << msAmericanFlatVol.count() << " milliseconds)" << std::endl;
    //


    std::cout << "Closed Form European Call Price: " << ssvi.normalizedPrice(x, T, true)*std::sqrt(F*K)*std::exp(-r*T) << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanCallLocalVol(S,K,T,true,false,driftFun,localVolFun,discountRateFun,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanCallLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local Vol PDE European Call Price: " << europeanCallLocalVol.price() << " (" << msEuropeanCallLocalVol.count() << " milliseconds)" << std::endl;
    //

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanCallFlatVol(S,K,T,mu,r,sigma,true,false,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanCallFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE European Call Price: " << europeanCallFlatVol.price() << " (" << msEuropeanCallFlatVol.count() << " milliseconds)" << std::endl;
    //

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanCallLocalVol(S,K,T,true,true,driftFun,localVolFun,discountRateFun,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanCallLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local-Vol PDE American Call Price: " << americanCallLocalVol.price() << " (" << msAmericanCallLocalVol.count() << " milliseconds)" << std::endl;
    //

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanCallFlatVol(S,K,T,mu,r,sigma,true,true,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanCallFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE American Call Price: " << americanCallFlatVol.price() << " (" << msAmericanCallFlatVol.count() << " milliseconds)" << std::endl;
    //
    

}


int main() {

    //europeanFlatSigma(false);

    //europeanLocalVolatility(false);
    //europeanLocalVolatility(false);

    functionMdFile();



}