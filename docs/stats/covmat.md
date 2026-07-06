### Covariance matrix 

Implementation of the covariance/correlation matrix which wrap an Eigen's Matrix and have specific methods related to the object. For example, cholesky decomposition can directly be obtained from a valid correlation matrix object. A correlation matrix can also be obtained directly from the covariance matrix object. For more info, check the wiki page: https://en.wikipedia.org/wiki/Covariance_matrix


#### C++ snippet example

```cpp
int main() {

    Eigen::MatrixXd covMatrix(3,3);
    Eigen::MatrixXd corrMatrix(3,3);

    covMatrix << 0.0004026385246,-0.000005150756644,-0.0001141394503, 
                -0.000005150756644,0.00144691174,0.00003245725223, 
                -0.0001141394503,	0.00003245725223, 0.0005348029993;

    corrMatrix << 1.0, -0.006748268411, -0.245969, 
                -0.006748268411, 1.0, 0.03689723423, 
                -0.245969, 0.03689723423, 1.0;

    Eigen::VectorXd variances(3);  
    variances << 0.0004026385246,0.00144691174,0.0005348029993;

    nummcpp::stats::CovarianceMatrix cov1(covMatrix);
    nummcpp::stats::CovarianceMatrix cov2(CorrelationMatrix{corrMatrix},variances);

    std::cout << cov1.choleskyDecomposition() << std::end;
    //           1           0           0
    // -0.00674827    0.999977           0
    // -0.245969   0.0352382    0.968637
    std::cout << cov2.correlationMatrix().eigenvalues() << std::endl;
    //0.752183
    //0.998021
    //1.2498

}
```