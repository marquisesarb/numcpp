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

    double S = 100.0;   
    double K = 100.0;    
    double r = 0.05;     
    double q = 0.02;     
    double sigma = 0.20;
    double T = 2.0;     

    double callPrice = 13.5218;
    double putPrice = 7.9266;

    double F = S*std::exp((r-q)*T);
    double x = std::log(F/K); 
    double normalizedSigma = sigma*sqrt(T);
    std::cout << numcpp::cfin::blackScholesEuropeanNormalizedPrice(x,T,sigma, true)*std::exp(-r*T)*std::sqrt(F*K) << std::endl;
    // 13.5218
    std::cout << numcpp::cfin::blackScholesEuropeanNormalizedPrice(x,T,sigma, false)*std::exp(-r*T)*std::sqrt(F*K) << std::endl;
    // 7.9266
    std::cout << numcpp::cfin::blackScholesImpliedVolatility(modelPutPrice*std::exp(r*T)/std::sqrt(F*K),x,T,false) << std::endl;
    // 0.2
    std::cout << numcpp::cfin::blackScholesImpliedVolatility(modelCallPrice*std::exp(r*T)/std::sqrt(F*K),x,T, true)<< std::endl;
    // 0.2

}
```