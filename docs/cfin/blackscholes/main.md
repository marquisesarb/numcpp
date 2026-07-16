### Black-Sholes model

Black-Sholes Model refers to the formulation of a spot price $S_t$ following an geometric brownian motion which satisfies the following SDE: 

$$
dS_t/S_t = \mu_t dt + \sigma(S_t,t)dW_t
$$

Where $\mu_t$ can be extracted from the forward term structure of the spot price, as shown [here](driftderiv.md). C++ Implementation of the Back Scholes formula and its inversion for "implied volatility" are provided in [here](../../../include/numcpp/cfin/blackscholes/cf.hpp). For more details, see the following scientific papers: 

- [1] The pricing of options and corporate liabilities - Black and Scholes (1973) : https://www.jstor.org/stable/1831029
- [2] The pricing of commodity contract - Black (1976) : https://www.sciencedirect.com/science/article/abs/pii/0304405X76900246
- [3] Let's be rational - Jäckel (2015) : http://www.jaeckel.org/
- [4] Efficient analytic approximation of American Option Values - Barone-Adesi and Whaley (1987) : https://doi.org/10.1111/j.1540-6261.1987.tb02569.x
- [5] Pricing with a Smile, Bruno Dupire (1994) : https://scholar.google.com/citations?view_op=view_citation&hl=en&user=iXmfS1IAAAAJ&citation_for_view=iXmfS1IAAAAJ:W7OEmFMy1HYC

The Euler discretization (presented more in details [here]()) is also implemented [here](../../../include/numcpp/cfin/blackscholes/mc.hpp), note that the "Local Volatility" case is directly supported by the class. Furthermore, the engine to solve the model's PDE (presented more in details [here](pdesolver.md)) for different payoffs is implemented [here](../../../include/numcpp/cfin/blackscholes/pdesolver/interface.hpp). Note that the engine also supports the passing of a local volatility function. This can be done since the black-scholes model is a flat local volatility surface. 

#### C++ snippet example

##### Closed form formulas 

```cpp
int main() {

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
    // Closed Form European Put Price: 50.2866

    auto t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanLocalVol(S,K,T,false,false,driftFun,localVolFun,discountRateFun,N,M, false,10);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local Vol PDE European Put Price: " << europeanLocalVol.price() << " (" << msEuropeanLocalVol.count() << " milliseconds)" << std::endl;
    // Local Vol PDE European Put Price: 50.2268 (175 milliseconds)

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanFlatVol(S,K,T,mu,r,sigma,false,false,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE European Put Price: " << europeanFlatVol.price() << " (" << msEuropeanFlatVol.count() << " milliseconds)" << std::endl;
    // Flat-Vol PDE European Put Price: 50.2302 (57 milliseconds)

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanLocalVol(S,K,T,false,true,driftFun,localVolFun,discountRateFun,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local-Vol PDE American Put Price: " << americanLocalVol.price() << " (" << msAmericanLocalVol.count() << " milliseconds)" << std::endl;
    // Local-Vol PDE American Put Price: 51.4663 (204 milliseconds)

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanFlatVol(S,K,T,mu,r,sigma,false,true,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE American Put Price: " << americanFlatVol.price() << " (" << msAmericanFlatVol.count() << " milliseconds)" << std::endl;
    // Flat-Vol PDE American Put Price: 51.7251 (87 milliseconds)


    std::cout << "Closed Form European Call Price: " << ssvi.normalizedPrice(x, T, true)*std::sqrt(F*K)*std::exp(-r*T) << std::endl;
    // Closed Form European Call Price: 33.1184

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanCallLocalVol(S,K,T,true,false,driftFun,localVolFun,discountRateFun,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanCallLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local Vol PDE European Call Price: " << europeanCallLocalVol.price() << " (" << msEuropeanCallLocalVol.count() << " milliseconds)" << std::endl;
    // Local Vol PDE European Call Price: 32.9583 (176 milliseconds)

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla europeanCallFlatVol(S,K,T,mu,r,sigma,true,false,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msEuropeanCallFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE European Call Price: " << europeanCallFlatVol.price() << " (" << msEuropeanCallFlatVol.count() << " milliseconds)" << std::endl;
    // Flat-Vol PDE European Call Price: 32.9953 (57 milliseconds)

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanCallLocalVol(S,K,T,true,true,driftFun,localVolFun,discountRateFun,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanCallLocalVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Local-Vol PDE American Call Price: " << americanCallLocalVol.price() << " (" << msAmericanCallLocalVol.count() << " milliseconds)" << std::endl;
    // Local-Vol PDE American Call Price: 39.4808 (204 milliseconds)

    t1 = std::chrono::high_resolution_clock::now();
    numcpp::cfin::BlackScholesPDESolverVanilla americanCallFlatVol(S,K,T,mu,r,sigma,true,true,N,M, false,10);
    t2 = std::chrono::high_resolution_clock::now();
    auto msAmericanCallFlatVol = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Flat-Vol PDE American Call Price: " << americanCallFlatVol.price() << " (" << msAmericanCallFlatVol.count() << " milliseconds)" << std::endl;
    // Flat-Vol PDE American Call Price: 39.2436 (88 milliseconds)

    

}
```
