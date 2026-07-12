#include <cassert>
#include <numcpp/cfin/blackscholes/cf.hpp>

bool isClose(double a, double b, double tol = 1e-3) {return std::abs(a - b) <= tol;}

void initialTestBlackScholes() {

    double S = 100.0;   
    double K = 100.0;    
    double r = 0.05;     
    double q = 0.02;     
    double sigma = 0.20;
    double T = 2.0;     

    double callPrice = 13.5218;
    double putPrice = 7.9266;

    double F = S*std::exp((r-q)*T);
    double x = std::log(F/K); 
    double normalizedSigma = sigma*sqrt(T);
    double modelCallPrice = numcpp::cfin::blackScholesEuropeanNormalizedPrice(x,T,sigma, true)*std::exp(-r*T)*std::sqrt(F*K);
    double modelPutPrice = numcpp::cfin::blackScholesEuropeanNormalizedPrice(x,T,sigma, false)*std::exp(-r*T)*std::sqrt(F*K);
    double impliedVolPut = numcpp::cfin::blackScholesImpliedVolatility(modelPutPrice*std::exp(r*T)/std::sqrt(F*K),x,T,false);
    double impliedVolCall = numcpp::cfin::blackScholesImpliedVolatility(modelCallPrice*std::exp(r*T)/std::sqrt(F*K),x,T, true);
    assert(isClose(modelCallPrice, callPrice, 1e-4));
    assert(isClose(modelPutPrice, putPrice, 1e-4));
    assert(isClose(impliedVolPut, .2, 1e-3));
    assert(isClose(impliedVolCall, .2, 1e-3));

}

void testTemplateHighDimensionLetsBeRational(double sigmaMin, double sigmaMax, double xMin, double xMax) {

    double sigmaDelta = (sigmaMax - sigmaMin)/100.0; 
    double xDelta = (xMax - xMin)/100.0; 
    double sigma = sigmaMin, x = xMin;
    std::vector<double> sigmas, moneyness; 
    double priceCall, estSigmaCall,pricePut, estSigmaPut;
    for (int i = 0; i<100; i++) {
        sigmas.push_back(sigma);
        moneyness.push_back(x);
        x += xDelta; 
        sigma += sigmaDelta;
    }
    
    for (double s : sigmas) {

        for (double xx: moneyness) {

            priceCall = numcpp::cfin::letsberational::_getNormalizedPrice(xx,s,true);
            pricePut = numcpp::cfin::letsberational::_getNormalizedPrice(xx,s,false);

            if (priceCall>numcpp::cfin::letsberational::_getNormalizedIntrisicValue(xx, true)) {
                estSigmaCall = numcpp::cfin::letsberational::_getNormalizedVolatility(priceCall, xx, true);
                if (!std::isnan(estSigmaCall)) {
                    assert(isClose((estSigmaCall-s)/s,0.0,1e-1));
                }
            }

            if (pricePut>numcpp::cfin::letsberational::_getNormalizedIntrisicValue(xx, false)) {
                estSigmaPut = numcpp::cfin::letsberational::_getNormalizedVolatility(pricePut, xx, false);
                if (!std::isnan(estSigmaPut)) {
                    assert(isClose((estSigmaPut-s)/s,0.0,1e-1));
                }
            }         
        }
    }
}

void testLetsBeRationalMainTest() {

    testTemplateHighDimensionLetsBeRational(10e-7,1.22,0,3);
    testTemplateHighDimensionLetsBeRational(10e-5,7.07,0,16);
    testTemplateHighDimensionLetsBeRational(10e-5,0.18,0,10e-5);

}

int main() {

    initialTestBlackScholes();
    testLetsBeRationalMainTest();
    return 0;
}
