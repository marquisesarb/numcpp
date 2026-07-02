#pragma once 
#include <numcpp/stats/distributions/gamma.hpp>

namespace numcpp::stats {

    struct ChiSquaredDistribution final: public GammaDistribution {


        ChiSquaredDistribution(double nu): GammaDistribution(nu/2.0,2.0){}

        double nu() const {return 2.0*k;}

    };

}