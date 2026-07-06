### Bisection Algorithm 

See wiki : https://en.wikipedia.org/wiki/Bisection_method
or mathematiczal optimization from Indranil Ghosh : https://indrag49.github.io/Numerical-Optimization/

The algortihm uses the mid-point method to reduce the range [a,b] by half and chooses which bound to replace with respect to the sign of $f(a)f(x)$, where $x=\frac{a+b}{2}$. If $f(a)f(x)$ is negative, then $b$ can be replace by x and an other iteration can be done.

#### C++ snippet example

```cpp
double f(double x) { return x * x - 2.0; };
double df(double x) { return 2.0 * x; };

int main() {

    numcpp::optim::Bisection bi(0.0, 2.0, f); 

    std::cout << bi.x << std::endl;
    // 1.41421
    std::cout << bi.fx << std::endl;
    // 6.7546e-13
    std::cout << bi.iter << std::endl;
    // 39 
}
```