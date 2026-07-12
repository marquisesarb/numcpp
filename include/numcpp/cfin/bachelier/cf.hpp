#pragma once
#include <cmath> 
#include <numcpp/stats/distributions/normal.hpp>
 
namespace numcpp::cfin {

    inline double bachelierEuropeanPrice(double F, double K, double t, double sigma, bool isCall) {
        
        double putCallFlag = isCall ? 1.0 : -1.0; 
        double d = (F-K)/(sigma*std::sqrt(t));
        numcpp::stats::NormalDistribution normal{};
        double Nd = normal.cdf(d*putCallFlag);
        double nd = normal.pdf(d*putCallFlag);
        return putCallFlag * (F-K) * Nd + sigma*std::sqrt(t)*nd;
    }

    inline double bachelierImpliedVolatility(double p, double F, double K, double t, bool isCall) {

        double otherOption =
            isCall
            ? p - (F - K)   
            : p + (F - K); 

        double straddlePrice = p + otherOption;
        double payoffStraddleRatio = (F-K)/straddlePrice; 
        double eta = payoffStraddleRatio/std::atanh(payoffStraddleRatio);
        double sqrteta = std::sqrt(eta);

        const double PI = 3.14159265358979323846;

        std::array<double, 8> a = {
            3.994961687345134e-1, 
            2.100960795068497e+1, 
            4.980340217855084e+1, 
            5.988761102690991e+2, 
            1.848489695437094e+3,
            6.106322407867059e+3,
            2.493415285349361e+4,
            1.266458051348246e+4};

        std::array<double, 10> b = {
            1.0, 
            4.990534153589422e+1, 
            3.093573936743112e+1, 
            1.495105008310999e+3, 
            1.323614537899738e+3,
            1.598919697679745e+4,
            2.392008891720782e+4, 
            3.608817108375034e+3, 
            -2.067719486400926e+2, 
            1.174240599306013e+1};

        double hnum = 0.0, hdenom = 0.0; 

        for (size_t i = 0; i<10; i++) {

            double etapow = (i==0) ? 1.0 : (i==1) ? eta : (i==2) ? eta*eta : std::pow(eta, double(i));
            hnum += (i>7) ? 0.0 : a[i]*etapow;
            hdenom += b[i]*etapow;
        }

        return straddlePrice * std::sqrt(0.5*PI/t) * sqrteta*hnum/hdenom;
    }

}