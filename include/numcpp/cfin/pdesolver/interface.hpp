#pragma once 
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <numcpp/solvers/tridiagthomas.hpp>

namespace numcpp::cfin {

    class PDESolverEngineInterface {

        public: 
            Eigen::VectorXd logMoneynessVector; 
            Eigen::MatrixXd spotPriceMatrix; 
            Eigen::MatrixXd optionValueMatrix; 

            PDESolverEngineInterface(int N, int M) {

                logMoneynessVector = Eigen::VectorXd::Zero(2*M+1); 
                spotPriceMatrix = Eigen::MatrixXd::Zero(2*M+1, N);
                optionValueMatrix = Eigen::MatrixXd::Zero(2*M+1, N);

            }

            virtual ~PDESolverEngineInterface() = default;

            double price() const {

                return optionValueMatrix.col(0)(size_t((logMoneynessVector.size()-1)/2));
            }
            
        protected:

            void computeLogMoneynessVector(double S, double dx) {

                size_t M = (logMoneynessVector.size()-1)/2;
                logMoneynessVector = Eigen::VectorXd::Zero(2*M+1);
                size_t mid = M;
                logMoneynessVector(mid) = 0.0;
                for (size_t i = 1; i<=mid; i++) {
                    logMoneynessVector(mid-i) = logMoneynessVector(mid-i+1)+dx;
                    logMoneynessVector(mid+i) = logMoneynessVector(mid+i-1)-dx;

                }
            }

            void computeSpotPriceMatrix(double S, double dt, const std::function<double(double)>& driftFunction) {

                double t = 0;
                double F = S;
                size_t M = (spotPriceMatrix.rows()-1)/2;
                spotPriceMatrix(M,0) = S;

                for (size_t i = 1; i<=M; i++) {
  
                    spotPriceMatrix(M-i,0) = F*std::exp(logMoneynessVector(M-i));
                    spotPriceMatrix(M+i,0) = F*std::exp(logMoneynessVector(M+i));

                }
                
                for (size_t i = 1; i<spotPriceMatrix.cols(); i++) {

                    t += dt; 
                    F = F*std::exp(driftFunction(t)*dt);
                    spotPriceMatrix(M,i) = F;

                    for (size_t u = 1; u<=M; u++) {
  
                        spotPriceMatrix(M-u,i) = F*std::exp(logMoneynessVector(M-u));
                        spotPriceMatrix(M+u,i) = F*std::exp(logMoneynessVector(M+u));

                    }

                }
            }

            static Eigen::SparseMatrix<double> flatSigmaImplicitMatrix(double dx, double dt, double sigma, double mu, double r, int M) {

                //Eigen::SparseMatrix<double> A(2*M+1, 2*M+1);
                Eigen::SparseMatrix<double> A(2*M+1, 2*M+1);
                std::vector<Eigen::Triplet<double>> coefs;
                //double mid = 1.0+dt*r + sigma*sigma*dt/(dx*dx); 
                //double up = -.5*sigma*sigma*dt/(dx*dx) - .5*(mu - .5*sigma*sigma)*dt/dx; 
                //double down = -.5*sigma*sigma*dt/(dx*dx) + .5*(mu - .5*sigma*sigma)*dt/dx; 

                double mid = 1.0+dt*r + sigma*sigma*dt/(dx*dx); 
                double up = -.5*sigma*sigma*dt/(dx*dx) - .5*(mu -.5*sigma*sigma)*dt/dx; 
                double down = -.5*sigma*sigma*dt/(dx*dx) + .5*(mu - .5*sigma*sigma)*dt/dx; 

                coefs.push_back({Eigen::Triplet<double>(0,0,2*up)});
                coefs.push_back({Eigen::Triplet<double>(0,1,2*mid-up)});
                coefs.push_back({Eigen::Triplet<double>(0,2,2*down-mid)});
                coefs.push_back({Eigen::Triplet<double>(0,3,-down)});
                coefs.push_back({Eigen::Triplet<double>(2*M,2*M,1.0+dt*r)});

                for (size_t i = 1; i<2*M; i++) {

                    coefs.push_back({Eigen::Triplet<double>(i,i-1,up)});
                    coefs.push_back({Eigen::Triplet<double>(i,i,mid)});
                    coefs.push_back({Eigen::Triplet<double>(i,i+1,down)});
                }

                A.setFromTriplets(coefs.begin(), coefs.end());
                return A;
            }

            static std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> flatSigmaCrankNicolsonMatrixes(double dx, double dt, double sigma, double mu, double r, int M) {


                Eigen::SparseMatrix<double> eM(2*M+1, 2*M+1);
                Eigen::SparseMatrix<double> iM(2*M+1, 2*M+1);
                std::vector<Eigen::Triplet<double>> eCoefs;
                std::vector<Eigen::Triplet<double>> iCoefs;
                double df = 1.0/(1.0+dt*r);
                double p = .5*sigma*sigma*dt/(dx*dx);

                double emid = df*(1-p);
                double eup = df*(.5*p + .25*(mu - .5*sigma*sigma)*dt/dx); 
                double edown = df*(.5*p - .25*(mu - .5*sigma*sigma)*dt/dx); 

                //double imid = 1.0+dt*r + .5*sigma*sigma*dt/(dx*dx); 
                double imid = 1.0+.5*dt*r + .5*sigma*sigma*dt/(dx*dx); 
                double iup = -.25*sigma*sigma*dt/(dx*dx) - .25*(mu - .5*sigma*sigma)*dt/dx; 
                double idown = -.25*sigma*sigma*dt/(dx*dx) + .25*(mu - .5*sigma*sigma)*dt/dx; 

                eCoefs.push_back({Eigen::Triplet<double>(0,0,2*eup)});
                eCoefs.push_back({Eigen::Triplet<double>(0,1,2*emid-eup)});
                eCoefs.push_back({Eigen::Triplet<double>(0,2,2*edown-emid)});
                eCoefs.push_back({Eigen::Triplet<double>(0,3,-edown)});
                eCoefs.push_back({Eigen::Triplet<double>(2*M,2*M,df)});

                iCoefs.push_back({Eigen::Triplet<double>(0,0,2*iup)});
                iCoefs.push_back({Eigen::Triplet<double>(0,1,2*imid-iup)});
                iCoefs.push_back({Eigen::Triplet<double>(0,2,2*idown-imid)});
                iCoefs.push_back({Eigen::Triplet<double>(0,3,-idown)});
                iCoefs.push_back({Eigen::Triplet<double>(2*M,2*M,1.0/df)});

                for (size_t i = 1; i<2*M; i++) {

                    eCoefs.push_back({Eigen::Triplet<double>(i,i-1,eup)});
                    eCoefs.push_back({Eigen::Triplet<double>(i,i,emid)});
                    eCoefs.push_back({Eigen::Triplet<double>(i,i+1,edown)});

                    iCoefs.push_back({Eigen::Triplet<double>(i,i-1,iup)});
                    iCoefs.push_back({Eigen::Triplet<double>(i,i,imid)});
                    iCoefs.push_back({Eigen::Triplet<double>(i,i+1,idown)});
                }

                eM.setFromTriplets(eCoefs.begin(), eCoefs.end());
                iM.setFromTriplets(iCoefs.begin(), iCoefs.end());
                return std::make_pair(eM,iM);
            }

            Eigen::SparseMatrix<double> localImplicitMatrix(
                double T, double dx, double dt, int M,
                const std::function<double(double)>& driftFunction, 
                const std::function<double(double,double)>& localVolatilityFunction, 
                const std::function<double(double)>& discountRateFunction) const {

                Eigen::SparseMatrix<double> A(2*M+1, 2*M+1);
                std::vector<Eigen::Triplet<double>> coefs;
                double r = discountRateFunction(T); 
                double sigma = localVolatilityFunction(-logMoneynessVector(0),T);
                double mu = driftFunction(T);
                double mid = 1.0+dt*r + sigma*sigma*dt/(dx*dx); 
                double up = -.5*sigma*sigma*dt/(dx*dx) - .5*(mu - .5*sigma*sigma)*dt/dx; 
                double down = -.5*sigma*sigma*dt/(dx*dx) + .5*(mu - .5*sigma*sigma)*dt/dx; 

                coefs.push_back({Eigen::Triplet<double>(0,0,2*up)});
                coefs.push_back({Eigen::Triplet<double>(0,1,2*mid-up)});
                coefs.push_back({Eigen::Triplet<double>(0,2,2*down-mid)});
                coefs.push_back({Eigen::Triplet<double>(0,3,-down)});
                coefs.push_back({Eigen::Triplet<double>(2*M,2*M,1.0+dt*r)});

                for (size_t i = 1; i<2*M; i++) {

                    sigma = localVolatilityFunction(-logMoneynessVector(i),T);
                    mid = 1.0+dt*r + sigma*sigma*dt/(dx*dx); 
                    up = -.5*sigma*sigma*dt/(dx*dx) - .5*(mu - .5*sigma*sigma)*dt/dx; 
                    down = -.5*sigma*sigma*dt/(dx*dx) + .5*(mu - .5*sigma*sigma)*dt/dx; 
                    coefs.push_back({Eigen::Triplet<double>(i,i-1,up)});
                    coefs.push_back({Eigen::Triplet<double>(i,i,mid)});
                    coefs.push_back({Eigen::Triplet<double>(i,i+1,down)});
                }

                A.setFromTriplets(coefs.begin(), coefs.end());
                return A;
            }

            std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> localCrankNicolsonMatrixes(
                double T, double dx, double dt, int M,
                const std::function<double(double)>& driftFunction, 
                const std::function<double(double,double)>& localVolatilityFunction, 
                const std::function<double(double)>& discountRateFunction) const {


                double r = discountRateFunction(T); 
                double sigma = localVolatilityFunction(-logMoneynessVector(0),T);
                double mu = driftFunction(T);
                Eigen::SparseMatrix<double> eM(2*M+1, 2*M+1);
                Eigen::SparseMatrix<double> iM(2*M+1, 2*M+1);
                std::vector<Eigen::Triplet<double>> eCoefs;
                std::vector<Eigen::Triplet<double>> iCoefs;
                double df = 1.0/(1.0+dt*r);
                double p = .5*sigma*sigma*dt/(dx*dx);

                double emid = df*(1-p);
                double eup = df*(.5*p + .25*(mu - .5*sigma*sigma)*dt/dx); 
                double edown = df*(.5*p - .25*(mu - .5*sigma*sigma)*dt/dx); 

                //double imid = 1.0+dt*r + .5*sigma*sigma*dt/(dx*dx); 
                double imid = 1.0+.5*dt*r + .5*sigma*sigma*dt/(dx*dx); 
                double iup = -.25*sigma*sigma*dt/(dx*dx) - .25*(mu - .5*sigma*sigma)*dt/dx; 
                double idown = -.25*sigma*sigma*dt/(dx*dx) + .25*(mu - .5*sigma*sigma)*dt/dx; 

                eCoefs.push_back({Eigen::Triplet<double>(0,0,2*eup)});
                eCoefs.push_back({Eigen::Triplet<double>(0,1,2*emid-eup)});
                eCoefs.push_back({Eigen::Triplet<double>(0,2,2*edown-emid)});
                eCoefs.push_back({Eigen::Triplet<double>(0,3,-edown)});
                eCoefs.push_back({Eigen::Triplet<double>(2*M,2*M,df)});

                iCoefs.push_back({Eigen::Triplet<double>(0,0,2*iup)});
                iCoefs.push_back({Eigen::Triplet<double>(0,1,2*imid-iup)});
                iCoefs.push_back({Eigen::Triplet<double>(0,2,2*idown-imid)});
                iCoefs.push_back({Eigen::Triplet<double>(0,3,-idown)});
                iCoefs.push_back({Eigen::Triplet<double>(2*M,2*M,1.0/df)});

                for (size_t i = 1; i<2*M; i++) {

                    sigma = localVolatilityFunction(-logMoneynessVector(i),T);
                    imid = 1.0+.5*dt*r + .5*sigma*sigma*dt/(dx*dx); 
                    iup = -.25*sigma*sigma*dt/(dx*dx) - .25*(mu - .5*sigma*sigma)*dt/dx; 
                    idown = -.25*sigma*sigma*dt/(dx*dx) + .25*(mu - .5*sigma*sigma)*dt/dx; 
                    emid = df*(1-p);
                    eup = df*(.5*p + .25*(mu - .5*sigma*sigma)*dt/dx); 
                    edown = df*(.5*p - .25*(mu - .5*sigma*sigma)*dt/dx); 

                    eCoefs.push_back({Eigen::Triplet<double>(i,i-1,eup)});
                    eCoefs.push_back({Eigen::Triplet<double>(i,i,emid)});
                    eCoefs.push_back({Eigen::Triplet<double>(i,i+1,edown)});

                    iCoefs.push_back({Eigen::Triplet<double>(i,i-1,iup)});
                    iCoefs.push_back({Eigen::Triplet<double>(i,i,imid)});
                    iCoefs.push_back({Eigen::Triplet<double>(i,i+1,idown)});
                }

                eM.setFromTriplets(eCoefs.begin(), eCoefs.end());
                iM.setFromTriplets(iCoefs.begin(), iCoefs.end());
                return std::make_pair(eM,iM);
            }



    };
}