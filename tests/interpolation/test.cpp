#include <cassert>
#include <numcpp/interpolation/linear.hpp>
#include <numcpp/interpolation/cubicspline.hpp>

int main() {

    std::map<double, double> testData = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 5.0},
        {4.0, 7.0},
        {5.0, 11.0}
    };
    
    numcpp::interpolation::LinearInterpolation linear(testData);
    
    assert(linear.evaluate(1.0) == 2.0);
    assert(linear.evaluate(3.0) == 5.0);
    assert(linear.evaluate(5.0) == 11.0);

    
    numcpp::interpolation::CubicSpline cubic(testData);
    
    assert(cubic.evaluate(1.0) == 2.0);
    assert(cubic.evaluate(3.0) == 5.0);
    assert(cubic.evaluate(5.0) == 11.0);

    assert(cubic.evaluateSecondDerivative(1.0) == 0.0);
    assert(cubic.evaluateSecondDerivative(5.0) == 0.0);



    return 0;
}

