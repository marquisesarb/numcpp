#pragma once 
#include "Eigen/Core"
#include <Eigen/Dense>
#include <numcpp/stats/distributions/tstudent.hpp>

namespace numcpp::reg {

    class OLS {

        public: 
            bool useIntercept; 
            Eigen::VectorXd betas; 
            Eigen::VectorXd residuals; 
            Eigen::MatrixXd invXtX;
            double sampleVariance;
            double rSquared;

            OLS(const Eigen::VectorXd& Y,const Eigen::MatrixXd& X, bool fitIntercept) {

                useIntercept = fitIntercept;
                Eigen::MatrixXd X_ = addIntercept(X, fitIntercept); 
                invXtX = inverseXtX(X_);
                betas = invXtX*X_.transpose()*Y;
                residuals = Y - X_*betas;
                double rss = residuals.dot(residuals);
                sampleVariance = rss / double(Y.size() - betas.size());
                double ymean = Y.mean();
                double tss = 0.0;
                for (size_t i = 0; i < Y.size(); ++i){double d = Y[i] - ymean;tss += d * d;}
                rSquared = 1.0 - rss / tss;
            }
            ~OLS() = default;

            int numberParameters() const {return betas.size();}

            int numberObservations() const {return residuals.size();}

            int degreesFreedom() const {return numberObservations()-numberParameters();}

            Eigen::VectorXd standardErrors() const {return (sampleVariance * invXtX.diagonal()).array().sqrt();}

            Eigen::VectorXd confidenceInterval(double quantile, bool upper) {

                double tcrit = numcpp::stats::TStudentDistribution{double(degreesFreedom())}.invCdf(0.5 + quantile / 2.0);
                //double tcrit = probability::invCdf(probability::Student{double(degreesFreedom())},0.5 + quantile / 2.0);
                Eigen::VectorXd ci(numberParameters());
                Eigen::VectorXd se = standardErrors();
                for(int j = 0; j < numberParameters(); ++j) {
                    ci(j) = betas(j) + (upper ? tcrit * se(j) : -tcrit * se(j));
                }

                return ci;
            }

            Eigen::VectorXd tstats() {

                Eigen::VectorXd ts(numberParameters()); 
                Eigen::VectorXd se = standardErrors();
                for(int j = 0; j < numberParameters(); ++j) {
                    ts(j) = betas(j)/se(j);
                }
                return ts;
            }

            double adjustedRSquared() const {return 1.0 - (1.0 - rSquared) * (numberObservations() - 1.0) / double(degreesFreedom());}

            Eigen::VectorXd pValues() {

                Eigen::VectorXd tstats_ = tstats();
                Eigen::VectorXd p(tstats_.size());

                for(int i = 0; i < tstats_.size(); ++i) {
                    p(i) = 2.0 * (1.0 - numcpp::stats::TStudentDistribution{double(degreesFreedom())}.cdf(std::abs(tstats_(i))));
                }

                return p;
            }

        private: 
      
            Eigen::MatrixXd addIntercept(const Eigen::MatrixXd& X, bool fitIntercept) {

                if (!fitIntercept) return X;
                Eigen::MatrixXd XX(X.rows(), X.cols() + 1);
                XX.col(0) = Eigen::VectorXd::Ones(X.rows());
                XX.block(0, 1, X.rows(), X.cols()) = X;
                return XX;
            }

            Eigen::MatrixXd inverseXtX(const Eigen::MatrixXd& X) {

                Eigen::MatrixXd XtX = X.transpose() * X;
                Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
                return ldlt.solve(Eigen::MatrixXd::Identity(XtX.rows(), XtX.rows()));
            }

            Eigen::MatrixXd projectionMatrix(const Eigen::MatrixXd& X) {

                Eigen::MatrixXd XtX = X.transpose() * X;
                Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
                return X * ldlt.solve(X.transpose());
            }


    };

}