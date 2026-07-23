#include <cassert>
#include <numcpp/interpolation/cubicspline.hpp>
#include <iostream>

bool isClose(double a, double b, double eps) { return std::abs(a-b)<eps;}


int main() {

    std::map<double, double> testData = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 5.0},
        {4.0, 7.0},
        {5.0, 11.0}
    };

    numcpp::interpolation::CubicSpline l1(testData); 

    numcpp::interpolation::CubicSpline l2 = numcpp::interpolation::CubicSpline(); 

    for (const auto& [k,v] : testData) {l2.insert(k, v);}

    assert(l1.size()==5);
    assert(l1.size()==l2.size()); 

    assert(std::isnan(l1.evaluate(.5)));
    assert(std::isnan(l1.evaluate(5.5)));
    assert(std::isnan(l2.evaluate(.5)));
    assert(std::isnan(l2.evaluate(5.5)));

    for (const auto& [k,v]: testData) {

        assert(isClose(l1.evaluate(k),testData.at(k), 1e-10));
        assert(isClose(l2.evaluate(k),testData.at(k), 1e-10));
    }

    std::cout << l1.evaluateSecondDerivative(1.0) << std::endl;
    std::cout << l1.evaluateSecondDerivative(5.0) << std::endl;
    assert(isClose(l1.evaluateSecondDerivative(1.0),0.0, 1e-10));
    assert(isClose(l1.evaluateSecondDerivative(5.0),0.0, 1e-10));
    assert(isClose(l2.evaluateSecondDerivative(1.0),0.0, 1e-10));
    assert(isClose(l2.evaluateSecondDerivative(5.0),0.0, 1e-10));

    assert(l1.contains(2.0)==true);
    l1.erase(2.0); 
    assert(l1.contains(2.0)==false);
    assert(l1.size()==4);
    assert(l1.index(2.0)==-1); 

    for (const auto& [k,v]: testData) {

        if (k!=2.0) {
            assert(isClose(l1.evaluate(k),testData.at(k), 1e-10));
        } 
    }

    assert(isClose(l1.evaluateSecondDerivative(1.0),0.0, 1e-10));
    assert(isClose(l1.evaluateSecondDerivative(5.0),0.0, 1e-10));

    l1.insert(2.5, 3.5); 
    assert(l1.size()==5);

    for (const auto& [k,v]: testData) {

        if (k!=2.0) {

            assert(isClose(l1.evaluate(k),testData.at(k), 1e-10));
        } 
    }

    assert(isClose(l1.evaluateSecondDerivative(1.0),0.0, 1e-10));
    assert(isClose(l1.evaluateSecondDerivative(5.0),0.0, 1e-10));

    assert(l1.evaluate(2.5)==3.5);
}