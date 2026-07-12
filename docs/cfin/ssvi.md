### Surface Stochastic Volatility Inspired 

Implementation of the Surface SVI (Power-Law) as defined from the following scientific paper: 

- [1] Arbitrage-free SVI volatility surfaces, Jim Gatheral and Antoine Jacquier : https://arxiv.org/abs/1204.0646

The Power-Law SSVI formalutation of the total variance reads: 

$$
w(k,t) = \frac{\theta_t}{2}\left(1 + \rho \varphi_t k + \sqrt{(\varphi_t k + \rho)^2 + (1-\rho^2)}\right)
$$

with $\rho \in [0,1]$ and $\theta_t : \mathbb{R}_+ \rightarrow \mathbb{R}_+$ represents the term structure of the ATM total variance. Also the Power-Law function is given by: 

$$
\varphi_t = \eta \theta_t^{-\gamma}
$$

with $\eta > 0$ and $\gamma \in (0,1)$. The total variance skew and ATM skew are given by: 

$$
\frac{\partial w}{\partial k} = \frac{\theta_t}{2}\left( \rho\varphi_t + \frac{\varphi_t (\varphi_t k+\rho)}{\sqrt{(\varphi_t k + \rho)^2 + (1-\rho^2)}} \right) \Longleftrightarrow \frac{\partial w}{\partial k}\mid_{k=0} = \theta_t \varphi_t \rho
$$

Also, the second derivative with respect to log-moneyness is given by: 

$$
\frac{\partial^2 w}{\partial k^2} = \frac{\theta_t\varphi_t^2(1-\rho^2)}{2\left(\sqrt{(\varphi_t k + \rho)^2 + (1-\rho^2)}\right)^3}
$$


#### C++ snippet example

```cpp
int main() {

    numcpp::cfin::SSVIPowerLawFlatATM ssvi(-0.5,0.46,0.77,1.0); 
    std::cout << glq.roots.transpose() << std::endl;
    // 0.137793 0.729455  1.80834  3.40143   5.5525  8.33015  11.8438  16.2793  21.9966  29.9207
    std::cout << glq.weights.transpose() << std::endl;
    // 0.308441     0.40112    0.218068   0.0620875  0.00950152 0.000753008 2.82592e-05 4.24931e-07 1.83956e-09 9.91183e-13

    std::cout << ssvi.totalVariance(2.0,0.1) << std::endl;
    // 0.164217944
    std::cout << ssvi.totalVarianceSkew(2.0,0.1) << std::endl;
    // 0.06566428559
    std::cout << ssvi.totaVarianceSecondDerivativeK(2.0,0.1) << std::endl;
    // 0.002210303315

}
```

