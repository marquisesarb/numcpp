#include <functional>
#include <limits>
#include <numcpp/cfin/blackscholes/baw.hpp>
#include <numcpp/stats/distributions/normal.hpp>

namespace numcpp::cfin::bawtools {

    double getQ(double mu, double sigma, double r, double t, bool isCall) {

        double putCallFlag = isCall ? 1.0 : -1.0; 
        double N = 2*mu/(sigma*sigma);
        double M = 2*r/(sigma*sigma);
        return .5*(1-N+putCallFlag*std::sqrt((N-1)*(N-1)+4*M/(1-getDiscountFactor(t,r))));
    }

    double getQInf(double mu, double sigma, double r, bool isCall) {

        double putCallFlag = isCall ? 1.0 : -1.0; 
        double N = 2*mu/(sigma*sigma);
        double M = 2*r/(sigma*sigma);
        return .5*(1-N+putCallFlag*std::sqrt((N-1)*(N-1)+4*M));
    }

    double getDiscountFactor(double t, double r) {

        return std::exp(-r*t);
    }

    double _getExercisePremium(double F, double mu, double K, double t, double sigma, double r, bool isCall, double E) {

        if (mu >= r and isCall){return 0.0;}
        double putCallFlag = isCall ? 1.0 : -1.0;
        double qq_ = getQ(mu,sigma,r,t,isCall);
        double delta = europeanDelta(E*std::exp(mu*t),mu,K,t,sigma,r,isCall);
        return putCallFlag * E * std::pow(F / E, qq_) * (1 - putCallFlag * delta) / qq_;
    }

    double _getPrice(double F, double mu, double K, double t, double sigma, double r, bool isCall, double E) {

        double putCallFlag = isCall ? 1.0 : -1.0;
        double euroPrice = europeanPrice(F*std::exp(mu*t),K,t,sigma,r,isCall);

        if (mu >= r and isCall){ 

            return euroPrice; 

        } else if (putCallFlag*F >= putCallFlag*E) { 

            return putCallFlag*(F-K); 

        } else {

            return euroPrice + _getExercisePremium(F,mu,K,t,sigma,r,isCall,E);

        }
    }

    double _getVega(double F, double mu, double K, double t, double sigma, double r, bool isCall, double E1, double E2, double eps) {

        double putCallFlag = isCall ? 1.0 : -1.0;
        double euroVega = europeanVega(F*std::exp(mu*t),K,t,sigma,r,isCall);

        if (mu >= r and isCall){ 

            return euroVega; 

        } else if (putCallFlag*F >= putCallFlag*E1) {

            return 0.0; 
        
        } else {

            double premium1 = _getExercisePremium(F,mu,K,t,sigma,r,isCall,E1);
            double premium2 = _getExercisePremium(F,mu,K,t,sigma+eps,r,isCall,E2);
            return euroVega + std::max((premium2-premium1)/eps,0.0);
        }
        
    }

    double _getImpliedVolatilityNewtonRaphson2(double p, double F, double K, double t, double r, bool isCall) {

        double x = numcpp::cfin::blackScholesImpliedVolatility(
            p/getDiscountFactor(t,r)/std::sqrt(F*K)
            , std::log(F/K)
            , t
            , isCall
        );
        double exercisePrice1, exercisePrice2, vega, calculatedPrice, epsilon = 0.001, fX, dfX, xStep = 0.0, threshold = 1e-5;
        int numberIterations_ = 0, maxIterations = 10;
        for (int i = 1; i <= maxIterations; ++i) {

            numberIterations_ += 1;
            exercisePrice1 = baroneAdesiWhaleyOptimalExercisePrice(0.0,K,t,x,r,isCall);
            exercisePrice2 = baroneAdesiWhaleyOptimalExercisePrice(0.0,K,t,x+epsilon,r,isCall);
            calculatedPrice = _getPrice(F,0.0,K,t,x,r,isCall,exercisePrice1); 
            vega = _getVega(F,0.0,K,t,x,r,isCall,exercisePrice1,exercisePrice2,epsilon);
            if (vega < __DBL_MIN__) { x *= 3; maxIterations += 10;  continue; }
            fX = std::log(1/calculatedPrice) - std::log(1/p);
            dfX = -vega/calculatedPrice;

            if (std::abs(fX) < threshold) break;
            
            if (std::abs(dfX) < std::sqrt(__DBL_MIN__)) {

                return std::numeric_limits<double>::quiet_NaN();
            }
            
            double xNew = x - fX / dfX;
            xStep = std::abs(xNew - x);
            x = xNew;

            if (xStep < threshold) break;

            if (i==maxIterations) break;
            
        }
        return x;
    }

    double _getImpliedVolatilityNewtonRaphson(double p, double F, double K, double t, double r, bool isCall) {

        double x0 = numcpp::cfin::blackScholesImpliedVolatility(
            p/getDiscountFactor(t,r)/std::sqrt(F*K)
            , std::log(F/K)
            , t
            , isCall
        );

        std::function<std::pair<double,double>(double)> target = [p,F,K,t,r,isCall](double x) {

            double exercisePrice1 = baroneAdesiWhaleyOptimalExercisePrice(0.0,K,t,x,r,isCall);
            double exercisePrice2 = baroneAdesiWhaleyOptimalExercisePrice(0.0,K,t,x+0.001,r,isCall);
            double calculatedPrice = _getPrice(F,0.0,K,t,x,r,isCall,exercisePrice1); 
            double vega = _getVega(F,0.0,K,t,x,r,isCall,exercisePrice1,exercisePrice2,0.001);
            return std::make_pair(std::log(1/calculatedPrice) - std::log(1/p), -vega/calculatedPrice);
        };

        numcpp::optim::NewtonRaphson res(x0, target,20, 1e-12);

        return res.x;
    }

    std::function<double(double)> _getTarget(double mu, double K, double t, double sigma, double r, bool isCall) {

        return [mu, K, t, sigma, r, isCall](double F) {

            double putCallFlag = isCall ? 1.0 : -1.0;
            double q = getQ(mu,sigma,r,t,isCall);
            double price = europeanPrice(F*std::exp(mu*t),K,t,sigma,r,isCall);
            double delta = europeanDelta(F*std::exp(mu*t),mu,K,t,sigma,r,isCall);
            return putCallFlag*(F-K) - price - putCallFlag*F*(1-putCallFlag*delta)/q;

        };
    }

    std::function<double(double)> _getTargetFirstDerivative(double mu, double K, double t, double sigma, double r, bool isCall) {

        return [mu, K, t, sigma, r, isCall](double F) {

            double putCallFlag = isCall ? 1.0 : -1.0;
            double q = getQ(mu,sigma,r,t,isCall);
            double delta = europeanDelta(F*std::exp(mu*t),mu,K,t,sigma,r,isCall);
            double gamma = europeanGamma(F*std::exp(mu*t),mu,K,t,sigma,r,isCall);
            return putCallFlag - delta*(1-1/q) - (putCallFlag-gamma*F)/q;

        };
    }

    double europeanPrice(double F, double K, double t, double sigma, double r, bool isCall) {
        return getDiscountFactor(t,r)*numcpp::cfin::blackScholesEuropeanNormalizedPrice(std::log(F/K),t,sigma,isCall)*std::sqrt(F*K);
    }

    double europeanDelta(double F, double mu, double K, double t, double sigma, double r, bool isCall) {

        int putCallFlag = isCall ? 1 : -1; 
        double d1 = (std::log(F / K) + 0.5 * sigma * sigma * t) / (sigma * std::sqrt(t)); 
        return getDiscountFactor(t,r)*putCallFlag*exp(mu*t)*numcpp::stats::NormalDistribution{}.cdf(putCallFlag*d1);
    }

    double europeanGamma(double F, double mu, double K, double t, double sigma, double r, bool isCall) {

        double d1 = (std::log(F / K) + 0.5 * sigma * sigma * t) / (sigma * std::sqrt(t)); 
        double driftSq = exp(mu*t)*std::exp(mu*t);
        return getDiscountFactor(t,r)*driftSq*numcpp::stats::NormalDistribution{}.pdf(d1)/(F*sigma*std::sqrt(t));
    }

    double europeanVega(double F, double K, double t, double sigma, double r, bool isCall) {

        double d1 = (std::log(F / K) + 0.5 * sigma * sigma * t) / (sigma * std::sqrt(t)); 
        return F*getDiscountFactor(t,r)*numcpp::stats::NormalDistribution{}.pdf(d1)*std::sqrt(t);
    }
}