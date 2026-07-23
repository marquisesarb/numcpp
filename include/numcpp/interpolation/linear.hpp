#pragma once 
#include <numcpp/interpolation/interface.hpp>

namespace numcpp::interpolation {

    class LinearInterpolation final: public Interpolation {
        public: 
            LinearInterpolation(const std::map<double, double>& data): Interpolation(data){}
            LinearInterpolation(){}
            ~LinearInterpolation() = default;
        
        protected: 
            double _evaluate(double x) override {
                int i = findIndex(x);
                std::pair<double,double> prev = (*this)[i-1];
                std::pair<double,double> next = (*this)[i];
                double x0 = prev.first, x1 = next.first;
                double y0 = prev.second, y1 = next.second;
                double t = (x - x0) / (x1 - x0);
                return y0 + t * (y1 - y0);
            }

            double _evaluateFirstDerivative(double x) override {

                int i = findIndex(x);
                std::pair<double,double> prev = (*this)[i-1];
                std::pair<double,double> next = (*this)[i];
                double x0 = prev.first, x1 = next.first;
                double y0 = prev.second, y1 = next.second;
                return (y1 - y0) / (x1 - x0);
            }

            double _evaluateSecondDerivative(double x) override { return 0.0;}
    };

}