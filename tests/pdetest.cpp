#include <iostream>
#include <numcpp/cfin/bspde/european.hpp>


int main() {


    numcpp::cfin::BlackScholesImplicitSchemeEuropeanEngine engine2(100,100,0.001-0.11,0.001,0.16,1.0,true,600,200);
    std::cout << engine2.price() << std::endl;

}