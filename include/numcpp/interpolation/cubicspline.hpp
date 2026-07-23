#pragma once 
#include <numcpp/interpolation/interface.hpp>

namespace numcpp::interpolation {


    class CubicSpline final : public Interpolation {

        public:
            CubicSpline(const std::map<double, double>& data): Interpolation(data){
                computedFlag = false;
                compute();
            }
            CubicSpline(){computedFlag = false;}
            ~CubicSpline() = default;

            void insert(std::pair<double,double> point) {

                Interpolation::insert(point);
                computedFlag = false;
            }

            void insert(double x, double y) {

                Interpolation::insert(x, y); 
                computedFlag = false;
            }
            
            void popBack() {

                Interpolation::popBack(); 
                computedFlag = false;
            }
            void erase(size_t i) {

                Interpolation::erase(i); 
                computedFlag = false;
            }

            void erase(double x) {

                Interpolation::erase(x); 
                computedFlag = false;
            }

        protected: 
            double _evaluate(double x) override {

                if (!computedFlag) compute();
                int i = findIndex(x);
                std::pair<double,double> prev = (*this)[i-1];
                double dx = x - prev.first;
                return prev.second + b_[i-1] * dx + c_[i-1] * dx * dx + d_[i-1] * dx * dx * dx;
            }

            double _evaluateFirstDerivative(double x) override {

                if (!computedFlag) compute();
                int i = findIndex(x);
                std::pair<double,double> prev = (*this)[i-1];
                double dx = x - prev.first;
                return b_[i-1] + 2.0*c_[i-1] * dx + 3.0*d_[i-1] * dx * dx;
            }

            double _evaluateSecondDerivative(double x) override {

                if (!computedFlag) compute();
                int i = findIndex(x);
                std::pair<double,double> prev = (*this)[i-1];
                double dx = x - prev.first;
                return 2.0*c_[i-1] + 6.0*d_[i-1] * dx;
            }

        private: 
            bool computedFlag;
            std::vector<double>  b_, c_, d_;

            void compute() {

                size_t n = xVector_.size() - 1; 
                std::vector<double> h(n);
                std::vector<double> alpha(n, 0.0);
                std::vector<double> l(n + 1, 0.0);
                std::vector<double> mu(n, 0.0);
                std::vector<double> z(n + 1, 0.0);

                l[0] = 1.0;  
                mu[0] = 0.0;
                z[0] = 0.0;

                for (size_t i = 0; i < n; ++i) {
                    h[i] = xVector_[i + 1] - xVector_[i];
                    if (i>=1) {
                        alpha[i] = (3.0 / h[i]) * (yVector_[i + 1] - yVector_[i]) - (3.0 / h[i - 1]) * (yVector_[i] - yVector_[i - 1]);
                        l[i] = 2.0 * (xVector_[i + 1] - xVector_[i - 1]) - h[i - 1] * mu[i - 1];
                        mu[i] = h[i] / l[i];
                        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
                    }
                }

                l[n] = 1.0;
                z[n] = 0.0;

                c_ = std::vector<double>(n+1, 0.0);
                b_ = std::vector<double>(n, 0.0);
                d_ = std::vector<double>(n, 0.0);
                
                for (int j = n - 1; j >= 0; --j) {
                    c_[j] = z[j] - mu[j] * c_[j + 1];
                    b_[j] = (yVector_[j + 1] - yVector_[j]) / h[j] - h[j] * (c_[j + 1] + 2.0 * c_[j]) / 3.0;
                    d_[j] = (c_[j + 1] - c_[j]) / (3.0 * h[j]);
                }

                computedFlag = true;
            }

    };


}