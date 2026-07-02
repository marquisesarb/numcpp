#include <numcpp/stats/distributions/uniform.hpp>
#include <numcpp/optim/brent.hpp>
#include <cassert>


std::pair<double,double> bruteMinimum(double a, double b, const std::function<double(double)>& f, int N) {

    double minVal=1e100, minX;

    for (int i = 0; i <= N; ++i) {
        double x = a + (b - a) * i / N;
        if (f(x)<minVal) {minVal = f(x); minX = x;}
    }
    return {minX, f(minX)};
}

std::function<double(double)> randomQuadratic(std::mt19937 & gen) {

    double a = numcpp::stats::UniformDistribution{}.sample(gen);

    auto f = [&](double x) {
        return (x - a)*(x - a);
    };

    return f;
}

void testBrentRandomQuadratic() {

    std::random_device rd;      
    std::mt19937 gen(rd());
    for (int i = 0; i<1000; i++) {
       
        std::function<double(double)> f = randomQuadratic(gen); 
        numcpp::optim::Brent res(f, -10.0, 10.0,50,1e-10); 
        
        assert(std::abs(res.fx) < 1e-6);
    }
}

void testOscilllaryFunction() {

    std::function<double(double)> f = [](double x) {
        return std::sin(x) + 0.1*x;
    };

    std::pair<double,double> brut = bruteMinimum(-5.0,5.0,f,10000);
    numcpp::optim::Brent brent(f,-5.0,5.0,50,1e-15);

    assert(std::abs(brut.first-brent.x)<1e-4);
    assert(std::abs(brut.second-brent.fx)<1e-5);
};

void testScaledQuadratic() {

    auto f = [](double x) {
        return 3.0*(x + 1.5)*(x + 1.5) + 5.0;
    };

   std::pair<double,double> brut = bruteMinimum(-5.0,5.0,f,10000);
    numcpp::optim::Brent brent(f,-5.0,5.0,50,1e-15);

    assert(std::abs(brut.first-brent.x)<1e-4);
    assert(std::abs(brut.second-brent.fx)<1e-5);
};

int main() {

    testBrentRandomQuadratic();
    testOscilllaryFunction();
    testScaledQuadratic();
}