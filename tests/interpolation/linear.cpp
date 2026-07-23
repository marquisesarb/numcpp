#include <cassert>
#include <numcpp/interpolation/linear.hpp>

int main() {

    std::map<double, double> testData = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 5.0},
        {4.0, 7.0},
        {5.0, 11.0}
    };

    numcpp::interpolation::LinearInterpolation l1(testData); 

    numcpp::interpolation::LinearInterpolation l2 = numcpp::interpolation::LinearInterpolation(); 

    for (const auto& [k,v] : testData) {l2.insert(k, v);}

    assert(l1.size()==5);
    assert(l1.size()==l2.size()); 

    double slope = (5.0-3.0)/(3.0-2.0);

    assert(l1.evaluate(5.0)==11.0);
    assert(l1.evaluate(3.0)==5.0);
    assert(l1.evaluate(5.0)==l2.evaluate(5.0)); 
    assert(l1.evaluate(3.0)==l2.evaluate(3.0)); 

    assert(std::isnan(l1.evaluate(.5)));
    assert(std::isnan(l1.evaluate(5.5)));
    assert(std::isnan(l2.evaluate(.5)));
    assert(std::isnan(l2.evaluate(5.5)));

    assert(l1.evaluate(2.5) == 3.0+slope*.5);
    assert(l1.evaluateFirstDerivative(2.5) == slope);
    assert(l1.evaluate(2.5)==l2.evaluate(2.5)); 
    assert(l1.evaluateFirstDerivative(2.5)==l2.evaluateFirstDerivative(2.5)); 

    assert(l1.contains(2.0)==true);
    l1.erase(2.0); 
    assert(l1.contains(2.0)==false);

    assert(l1.size()==4);
    assert(l1.index(2.0)==-1); 

    double slope2 = (5.0-2.0)/(3.0-1.0);
    assert(l1.evaluate(2.5) == 2.0+slope2*1.5);
    assert(l1.evaluateFirstDerivative(2.5) == slope2);


    l1.insert(2.5, 3.5); 
    double slope3 = (5.0-3.5)/(3.0-2.5);
    assert(std::abs(l1.evaluate(2.7) - (3.5+slope3*.2)) < 1e-10);
    assert(l1.evaluateFirstDerivative(2.7) == slope3);

}