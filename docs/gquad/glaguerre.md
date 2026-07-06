### Gauss-Laguerre 

See wiki : https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature

Note that weights and roots computations used Eigen's Adjoint solver object. 

#### C++ snippet example

```cpp
int main() {

    numcpp::gquad::GaussLaguerreQuadrature glq(10);
    std::cout << glq.roots.transpose() << std::endl;
    // 0.137793 0.729455  1.80834  3.40143   5.5525  8.33015  11.8438  16.2793  21.9966  29.9207
    std::cout << glq.weights.transpose() << std::endl;
    // 0.308441     0.40112    0.218068   0.0620875  0.00950152 0.000753008 2.82592e-05 4.24931e-07 1.83956e-09 9.91183e-13

}
```