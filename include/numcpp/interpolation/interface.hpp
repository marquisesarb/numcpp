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

            virtual ~Interpolation() = default;

            double evaluate(double x) const {

                if (x < getLowerBoundX() || x > getUpperBoundX()) return std::numeric_limits<double>::quiet_NaN();
                return _evaluate(x);
            }

            double evaluateFirstDerivative(double x) const{

                if (x < getLowerBoundX() || x > getUpperBoundX()) return std::numeric_limits<double>::quiet_NaN();
                return _evaluateFirstDerivative(x);
            }

            double evaluateSecondDerivative(double x) const{

                if (x < getLowerBoundX() || x > getUpperBoundX()) return std::numeric_limits<double>::quiet_NaN();
                return _evaluateSecondDerivative(x);
            }

            double getLowerBoundX() const { return xVector_.front();}
            double getUpperBoundX() const { return xVector_.back();}
            const std::vector<double>& getVectorX() const { return xVector_; }
            const std::vector<double>& getVectorY() const  { return yVector_; }

        protected: 
            virtual double _evaluate(double x) const = 0;
            virtual double _evaluateFirstDerivative(double x) const = 0; 
            virtual double _evaluateSecondDerivative(double x) const = 0; 

            int findIndex(double x) const {

                if (x==getLowerBoundX()) return 1; if (x==getUpperBoundX()) return xVector_.size();
                return std::distance(xVector_.begin(), std::lower_bound(xVector_.begin(), xVector_.end(), x));
            }
            double getX(int i) const { return xVector_[i];}
            double getY(int i) const { return yVector_[i];}

        private:
            std::vector<double> xVector_; 
            std::vector<double> yVector_; 

    };

}