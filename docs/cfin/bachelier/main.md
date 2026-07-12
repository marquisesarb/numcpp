### Bachelier model

Bachelier Model refers to the formulation of a spot price $S_t$ following an arithmetic brownian motion which satisfies the following SDE: 

$$
dS_t = \mu_t dt + \sigma(S_t,t)dW_t
$$

Where $\mu_t$ can be extracted from the forward term structure of the spot price, as shown [here](driftderiv.md). C++ Implementation of the Bachelier formula and its inversion for "implied volatility" are provided in [here](../../../include/numcpp/cfin/bachelier/cf.hpp). For more details, see the following scientific papers: 

- [1] Théorie de la spéculation - Bachelier (1900) : http://www.numdam.org/item?id=ASENS_1900_3_17__21_0
- [2] Numerical approximation of the implied volatility under arithmetic Brownian motion - Choi, Kim, Kwak (2007) : https://ssrn.com/abstract=990747

#### C++ snippet example

```cpp
int main() {

    double futurePrice = 100.0; 
    double sigma = 25;
    double timeToMaturity = 1.5;
    double strike = 115; 

    double callPrice = numcpp::cfin::bachelierEuropeanPrice(futurePrice, strike, timeToMaturity, sigma, true);
    double putPrice = numcpp::cfin::bachelierEuropeanPrice(futurePrice, strike, timeToMaturity, sigma, false);

    std::cout << callPrice << std::endl;
    // 6.15224
    std::cout << putPrice << std::endl;
    // 21.1522
    std::cout << numcpp::cfin::bachelierImpliedVolatility(callPrice, futurePrice, strike, timeToMaturity, true) << std::endl;
    // 25
    std::cout << numcpp::cfin::bachelierImpliedVolatility(putPrice, futurePrice, strike, timeToMaturity, false) << std::endl;
    // 25

}
```
