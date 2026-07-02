#include <numcpp/cfin/nss.hpp>
#include <cassert>

bool isClose(double a, double b, double eps) { return std::abs(a-b)<=eps;}

int main() {

    numcpp::cfin::NelsonSiegel ns{0.04,-0.018, 0.02,4.5}; 
    numcpp::cfin::Svensson sv{0.04,-0.018, 0.02,	-0.012,1,4.5};


    assert(isClose(ns.spotRate(.25),0.02302626974,1e-8));
    assert(isClose(ns.spotRate(6.75),0.03657322325,1e-8));

    assert(isClose(sv.spotRate(.25),0.02587233698,1e-8));
    assert(isClose(sv.spotRate(6.75),0.03673513498,1e-8));

    assert(isClose(ns.instantaneousForwardRate(.25),0.02402379564,1e-8));
    assert(isClose(ns.instantaneousForwardRate(6.75),0.04267756192,1e-8));

    assert(isClose(sv.instantaneousForwardRate(.25),0.02924495017,1e-8));
    assert(isClose(sv.instantaneousForwardRate(6.75),0.03612065003,1e-8));

}

