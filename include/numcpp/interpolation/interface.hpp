#pragma once 
#include <map>
#include <vector>
#include <limits>

namespace numcpp::interpolation {

    class Interpolation {
        public:
            Interpolation(const std::map<double, double>& data) {

                for(auto const& imap: data){xVector_.push_back(imap.first);yVector_.push_back(imap.second);}    
            }

            Interpolation() {}

            virtual ~Interpolation() = default;

            double evaluate(double x) {

                if (x < front().first || x > back().first) return std::numeric_limits<double>::quiet_NaN();
                return _evaluate(x);
            }

            double evaluateFirstDerivative(double x) {

                if (x < front().first || x > back().first) return std::numeric_limits<double>::quiet_NaN();
                return _evaluateFirstDerivative(x);
            }

            double evaluateSecondDerivative(double x) {

                if (x < front().first || x > back().first) return std::numeric_limits<double>::quiet_NaN();
                return _evaluateSecondDerivative(x);
            }

            std::pair<double,double> front() const { return std::make_pair(xVector_.front(), yVector_.front());}
            std::pair<double,double> back() const { return std::make_pair(xVector_.back(), yVector_.back());}
            std::pair<double,double> operator[](size_t i) const {return std::make_pair(xVector_[i], yVector_[i]);}
            bool empty() const {return xVector_.empty();}
            size_t size() const {return xVector_.size();}

            void insert(std::pair<double,double> point) {

                double x = point.first, y = point.second; 
                if (empty()) {xVector_.push_back(x); yVector_.push_back(y);}
                else {
                    std::pair<double,double> front_ = front(); 
                    std::pair<double,double> back_ = back(); 

                    if (x < front_.first) {xVector_.insert(xVector_.begin(),x); yVector_.insert(yVector_.begin(),y);}
                    else if (x > back_.first) {xVector_.push_back(x); yVector_.push_back(y);}
                    else {
                        auto it = std::lower_bound(xVector_.begin(), xVector_.end(), x);
                        size_t i = it - xVector_.begin();
                        if (x != xVector_[i]) {

                            xVector_.insert(xVector_.begin() + i, x);
                            yVector_.insert(yVector_.begin() + i, y);
                        }
                    }
                }
            }

            void insert(double x, double y) {insert(std::make_pair(x,y));}

            bool contains(double x) const {return (std::find(xVector_.begin(), xVector_.end(), x) != xVector_.end());}
            
            void popBack() {xVector_.pop_back();yVector_.pop_back();}
            void erase(size_t i) {xVector_.erase(xVector_.begin()+i); yVector_.erase(yVector_.begin()+i);}
            void erase(double x) {size_t i = index(x); if (i!=-1) erase(i);}

            size_t index(double x) const {

                auto it = std::lower_bound(xVector_.begin(), xVector_.end(), x);
                size_t i = it - xVector_.begin();
                return (xVector_[i] == x) ? i : -1;
            }

        protected: 
            virtual double _evaluate(double x) = 0;
            virtual double _evaluateFirstDerivative(double x) = 0; 
            virtual double _evaluateSecondDerivative(double x) = 0; 

            size_t findIndex(double x) const {

                if (x==front().first) return 1; if (x==back().first) return xVector_.size();
                return std::distance(xVector_.begin(), std::lower_bound(xVector_.begin(), xVector_.end(), x));
            }

            std::vector<double> xVector_; 
            std::vector<double> yVector_; 

    };

}