### Linear interpolation 

see wiki page : https://en.wikipedia.org/wiki/Linear_interpolation

#### C++ snippet example

```cpp
int main() {

    std::map<double, double> testData = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 5.0},
        {4.0, 7.0},
        {5.0, 11.0}
    };

    numcpp::interpolation::LinearInterpolation linear(testData);
    
    std::cout << linear.evaluate(1.5) << std::endl;
    // 2.5
    std::cout << linear.evaluateFirstDerivative(1.5) << std::endl;
    // 1
    std::cout << linear.evaluateSecondDerivative(3.5) << std::endl;
    // 0

}
```