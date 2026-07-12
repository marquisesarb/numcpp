#include <numcpp/cfin/bachelier/cf.hpp>
#include <cassert>

bool isClose(double a, double b, double tol = 1e-3) {return std::abs(a - b) <= tol;}

void testBachelierPrice() {

    double futurePrice = 100.0; 
    double sigma = 25;
    double timeToMaturity = 1.5;

    std::vector<double> strikes = {
        5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
        55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
        105, 110, 115, 120, 125, 130, 135, 140, 145, 150,
        155, 160, 165, 170, 175, 180, 185, 190, 195
    };

    std::vector<double> calls = {
        95.00810355, 90.01447127, 85.02524500, 80.04303039, 75.07168173,
        70.11673074, 65.18587699, 60.28950761, 55.44119414, 50.65809472,
        45.96117514, 41.37515770, 36.92811768, 32.65067482, 28.57477273,
        24.73209682, 21.15224297, 17.86080315, 14.87756934, 12.21506280,
        9.877569337, 7.860803154, 6.152242970, 4.732096821, 3.574772727,
        2.650674822, 1.928117685, 1.375157702, 0.9611751399, 0.6580947239,
        0.4411941367, 0.2895076115, 0.1858769938, 0.1167307369, 0.07168172529,
        0.0430303931, 0.02524499593, 0.0144712690, 0.008103551957
    };

    std::vector<double> puts = {
        0.008103551957, 0.0144712690, 0.02524499593, 0.0430303931, 0.07168172529,
        0.1167307369, 0.1858769938, 0.2895076115, 0.4411941367, 0.6580947239,
        0.9611751399, 1.375157702, 1.928117685, 2.650674822, 3.574772727,
        4.732096821, 6.152242970, 7.860803154, 9.877569337, 12.21506280,
        14.87756934, 17.86080315, 21.15224297, 24.73209682, 28.57477273,
        32.65067482, 36.92811768, 41.37515770, 45.96117514, 50.65809472,
        55.44119414, 60.28950761, 65.18587699, 70.11673074, 75.07168173,
        80.04303039, 85.02524500, 90.01447127, 95.00810355
    };

    for (size_t i = 0; i<strikes.size(); i++) {

        assert(isClose(numcpp::cfin::bachelierEuropeanPrice(futurePrice, strikes[i], timeToMaturity, sigma, true), calls[i],1e-8));
        assert(isClose(numcpp::cfin::bachelierEuropeanPrice(futurePrice, strikes[i], timeToMaturity, sigma, false), puts[i],1e-8));
    }

}

void testBachelierImpliedVol() {

    double futurePrice = 100.0;
    double timeToMaturity = 1.5;

    std::vector<double> strikes; 
    std::vector<double> sigmas;
    double strikeFactor = 0.05; 
    double sigmaFactor = 0.05; 
    for (int i = 0; i<41; i++) {
        strikes.push_back(strikeFactor*futurePrice); 
        sigmas.push_back(sigmaFactor*futurePrice); 
        strikeFactor += .05;
        sigmaFactor += .05;
    }

    for (double k: strikes) {

        for (double s: sigmas) {

            double call = numcpp::cfin::bachelierEuropeanPrice(futurePrice, k, timeToMaturity, s, true);
            double put = numcpp::cfin::bachelierEuropeanPrice(futurePrice, k, timeToMaturity, s, false);
            double ivcall = numcpp::cfin::bachelierImpliedVolatility(call, futurePrice, k, timeToMaturity, true);
            double ivput = numcpp::cfin::bachelierImpliedVolatility(put, futurePrice, k, timeToMaturity, false);
            

            if (put-std::max(-futurePrice+k,0.0)>1e-5 and call-std::max(futurePrice-k,0.0)>1e-5) {

                assert(isClose(ivcall, s, 1e-7));
                assert(isClose(ivput, s, 1e-7));
            }
            
        }
    }


}

int main() {

    testBachelierPrice(); 
    testBachelierImpliedVol();
    return 0;
}