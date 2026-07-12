#include <numcpp/cfin/ssvi.hpp>
#include <cassert>

bool isClose(double a, double b, double eps) { return std::abs(a-b)<=eps;}

int main() {

    numcpp::cfin::SSVIPowerLawFlatATM ssvi(-0.5,0.46,0.77,1.0); 

    assert(isClose(ssvi.totalVariance(-4.0,0.1), 0.8892560461, 1e-8)); 
    assert(isClose(ssvi.totalVarianceSkew(-4.0,0.1), -0.202757474, 1e-8)); 
    assert(isClose(ssvi.totaVarianceSecondDerivativeK(-4.0,0.1), 0.0001872933355, 1e-8)); 

    assert(isClose(ssvi.totalVariance(-2.0,0.1), 0.4844539759, 1e-8)); 
    assert(isClose(ssvi.totalVarianceSkew(-2.0,0.1), -0.2017235032, 1e-8)); 
    assert(isClose(ssvi.totaVarianceSecondDerivativeK(-2.0,0.1), 0.001286342193, 1e-8)); 

    assert(isClose(ssvi.totalVariance(.0,0.1), 0.1000, 1e-8)); 
    assert(isClose(ssvi.totalVarianceSkew(.0,0.1), -0.1354340407, 1e-8)); 
    assert(isClose(ssvi.totaVarianceSecondDerivativeK(.0,0.1), 0.2751356908, 1e-8)); 

    assert(isClose(ssvi.totalVariance(2.0,0.1), 0.164217944, 1e-8)); 
    assert(isClose(ssvi.totalVarianceSkew(2.0,0.1), 0.06566428559, 1e-8)); 
    assert(isClose(ssvi.totaVarianceSecondDerivativeK(2.0,0.1), 0.002210303315, 1e-8)); 

    assert(isClose(ssvi.totalVariance(4.0,0.1), 0.2976791797, 1e-8)); 
    assert(isClose(ssvi.totalVarianceSkew(4.0,0.1), 0.06724399819, 1e-8)); 
    assert(isClose(ssvi.totaVarianceSecondDerivativeK(4.0,0.1), 0.000246656128, 1e-8)); 

}