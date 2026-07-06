
// Implementation fo example here: https://www.youtube.com/watch?v=to82dv2SX28&time_continue=451&source_ve_path=MjM4NTE&embeds_referring_euri=https%3A%2F%2Fwww.google.com%2F

#include <iostream>
#include <numcpp/cfin/bspde/european.hpp>


int main() {

    //std::cout << blackScholesImplicitSchemeEuropeanPrice(100,100,0.05,0.01,0.2,2.0,true,400,10) << std::endl;
    numcpp::cfin::BlackScholesImplicitSchemeEuropeanEngine engine(100,100,0.001-0.11,0.001,0.16,1.0,true,252,79);
    std::cout << engine.price() << std::endl;

}