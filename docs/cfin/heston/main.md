### Heston model

Heston Model refers to the formulation of a spot price $S_t$ following an geometric brownian motion with CIR stochastic volatility which satisfies the following SDE: 

$$
dS_t/S_t = \mu_t dt + \sqrt{v_t}dW^S_t
$$

$$
dv_t = \kappa(\theta - v_t) + \eta \sqrt{v_t}dW^v_t
$$

with $<dW^S_t, dW^v_t> = \rho dt$

Where $\mu_t$ can be extracted from the forward term structure of the spot price, as shown [here](../blackscholes/driftderiv.md). C++ Implementation of the Heston characteristic function and its price formula (using Lewis framework and gaussian quadratures implemented [here](../../../include/numcpp/cfin/lewis.hpp)) are provided in [here](../../../include/numcpp/cfin/heston/cf.hpp). For more details, see the following scientific papers: 

- [1] The volatility surface, a practitioner’s guide - Gatheral (2006) : https://link.springer.com/article/10.1007/s11408-007-0072-4
- [2] The Heston Model and its Extensions in Matlab and C# - Fabrice Douglas Rouah (2013): https://onlinelibrary.wiley.com/doi/book/10.1002/9781118656471

#### C++ snippet example

##### Closed form formulas 

```cpp

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

    std::cout << heston.europeanNormalizedPriceGaussLaguerre(x,T,true,gl)*std::exp(-r*T)*sqrt(F*K) << std::endl;
    //13.2561

    std::cout << heston.impliedVolatilityGaussLaguerre(x,T,gl) << std::endl;
    //0.217483

    std::cout << heston.impliedDensityGaussLaguerre(x,T,gl) << std::endl;
    //1.33067

}
```

#### Quadratic exponential monte carlo

```cpp

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
    // Call MC price: 18.8929
    std::cout << "Call semi-closed form price: " << heston.europeanNormalizedPriceGaussLaguerre(x,T,true,gl)*sqrt(F*K) << std::endl;
    // Call semi-closed form price: 18.7514
    std::cout << "Put MC price: " << priceEuropeanOptionFromMatrix(sim.xMatrix, false, K) << std::endl;
    // Put MC price: 8.21682
    std::cout << "Put semi-closed form price: " << heston.europeanNormalizedPriceGaussLaguerre(x,T,false,gl)*sqrt(F*K) << std::endl;
    // Put semi-closed form price: 8.23433

    return 0;

}
```


