#pragma once 
#include <numcpp/interpolation/interface.hpp>

namespace numcpp::interpolation {

    class LinearInterpolation final: public Interpolation {
            public: 
                LinearInterpolation(const std::map<double, double>& data): Interpolation(data){}
                ~LinearInterpolation() = default;
            
            protected: 
                double _evaluate(double x) const override {
                    int i = findIndex(x);
                    double x0 = getX(i - 1), x1 = getX(i);
                    double y0 = getY(i - 1), y1 = getY(i);
                    double t = (x - x0) / (x1 - x0);
                    return y0 + t * (y1 - y0);
                }

                double _evaluateFirstDerivative(double x) const override {

                    int i = findIndex(x);
                    double x0 = getX(i - 1), x1 = getX(i);
                    double y0 = getY(i - 1), y1 = getY(i);
                    return (y1 - y0) / (x1 - x0);
                }

                double _evaluateSecondDerivative(double x) const override { return 0.0;}
        };

}