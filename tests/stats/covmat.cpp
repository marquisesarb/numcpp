#include <numcpp/stats/covmat.hpp>
#include <cassert>

bool isClose(double a, double b, double eps) { return std::abs(a-b)<eps;}

void assetIsMatrixesSimilar(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, double tol = 1e-9) {

    assert(A.rows() == B.rows() && A.cols() == B.cols());
    assert(A.isApprox(B, tol));
}

void testCorrelationMatrix3x3() {

    using namespace numcpp::stats; 

    Eigen::VectorXd corrCoeffs(3) ; 
    corrCoeffs << -0.006748268411, -0.2459694651, 0.03689723423; 

    Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Zero(3,3); 
    corrMatrix << 1.0, -0.006748268411,  -0.2459694651,
                -0.006748268411, 1.0, 0.03689723423,
                -0.2459694651, 0.03689723423, 1.0;
    
    Eigen::MatrixXd cholesky = Eigen::MatrixXd::Zero(3,3);  
    cholesky << 1.0, 0.0, 0.0,
                -0.006748268411, 1.0, 0.0,
                -0.2459694651, 0.03524, 0.9686;

    CorrelationMatrix corrMat(corrCoeffs); 

    assetIsMatrixesSimilar(corrMatrix,corrMat.matrix);

    assert(corrMat.matrix.rows()==3); 
    assert(corrMat.correlationCoefficient(0,0)==1); 
    assert(corrMat.correlationCoefficient(0,2)==-0.2459694651);

    assetIsMatrixesSimilar(cholesky,corrMat.choleskyDecomposition(), 1e-4);

}

void testCorrelationMatrix2x2() {

    using namespace numcpp::stats; 

    double corr = -0.006748268411;

    Eigen::MatrixXd corrMatrix(2,2); 
    corrMatrix << 1.0, -0.006748268411,
                -0.006748268411, 1.0;
    
    Eigen::MatrixXd cholesky(2,2); 
    cholesky << 1.0, 0.0, 
                -0.007, 1.0;

    CorrelationMatrix corrMat(corr); 
    assetIsMatrixesSimilar(corrMatrix,corrMat.matrix);

    assert(corrMat.matrix.rows()==2); 
    assert(corrMat.correlationCoefficient(0,0)==1); 
    assert(corrMat.correlationCoefficient(0,1)==-0.006748268411);
  
    assetIsMatrixesSimilar(cholesky,corrMat.choleskyDecomposition(), 1e-3);

}

void testCorrelationMatrix3x3MatrixConstructor() {

    using namespace numcpp::stats; 
    Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Zero(3,3); 
    corrMatrix << 1.0, -0.006748268411, -0.2459694651, 
        -0.006748268411, 1.0, 0.03689723423, 
        -0.2459694651, 0.03689723423, 1.0;

    CorrelationMatrix corrMat(corrMatrix);

    Eigen::MatrixXd cholesky = Eigen::MatrixXd::Zero(3,3); 
    cholesky <<  1.0, 0.0, 0.0, 
        -0.006748268411, 1.0, 0.0, 
        -0.2459694651, 0.03524, 0.9686; 
    

    assetIsMatrixesSimilar(corrMatrix,corrMat.matrix);

    assert(corrMat.correlationCoefficient(0,0)==1); 
    assert(corrMat.correlationCoefficient(0,2)==-0.2459694651);

    assetIsMatrixesSimilar(cholesky,corrMat.choleskyDecomposition(), 1e-4);


}

void testCorrelationMatrixNonSquareError() {

    using namespace numcpp::stats; 

    Eigen::MatrixXd corrMatrix(3,4); 
    corrMatrix << 1.0, -0.006748268411, -0.2459694651,1.0, 
        -0.006748268411, 1.0, 0.03689723423, 1.0, 
        -0.2459694651, 0.03689723423, 1.0, 1.0;

    assert(CorrelationMatrix{corrMatrix}.isValid()==false); 

    corrMatrix = Eigen::MatrixXd::Zero(3,3);
    corrMatrix << 1.0, -0.03689723423, -0.2459694651, 
        -0.006748268411, 1.0, 0.03689723423, 
        -0.2459694651, 0.03689723423, 1.0;

    assert(CorrelationMatrix{corrMatrix}.isValid()==false); 

    corrMatrix << 1.0, -0.03689723423, -0.2459694651, 
        -0.006748268411, 1.0, 0.03689723423, 
        -0.2459694651, 0.03689723423, 1.0;

    assert(CorrelationMatrix{corrMatrix}.isValid()==false); 
    
}

void testCorrelationMatrixNonSymmetricError() {

    using namespace numcpp::stats; 
    Eigen::MatrixXd corrMatrix(3,3);
    corrMatrix << 
        1.0, -0.03689723423, -0.2459694651, 
        -0.006748268411, 1.0, 0.03689723423, 
        -0.2459694651, 0.03689723423, 1.0;
    
    assert(CorrelationMatrix{corrMatrix}.isValid()==false); 

}

void testCorrelationMatrixNonPositiveSemiDefiniteError() {

    using namespace numcpp::stats; 
    Eigen::MatrixXd corrMatrix(3,3);
    corrMatrix << 1.0, -1.0, 0.0, 
        -1.0, 1.0, 0.0, 
        0.0, 0.0, 1.0;

    assert(CorrelationMatrix{corrMatrix}.isValid()==false); 

}

void testCorrelationMatrixInvalidCoefficients() {

    using namespace numcpp::stats; 
    Eigen::MatrixXd corrMatrix(2,2); 
    corrMatrix << 1.0, -2.0, 
        -2.0, 1.0;
    assert(CorrelationMatrix{corrMatrix}.isValid()==false); 

    assert(CorrelationMatrix::isNumberCoefficientsValid(4)==false); 

    assert(CorrelationMatrix::isNumberCoefficientsValid(6)==true); 
    
}

void testCovarianceMatrixConstructor() {

    using namespace numcpp::stats; 
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

    CovarianceMatrix cov1(covMatrix);
    assert(cov1.correlationMatrix().matrix.rows()==3);
    assert(cov1.matrix(0,1) == -0.000005150756644);
    assert(cov1.variance(1) == 0.00144691174);
    assetIsMatrixesSimilar(covMatrix,cov1.matrix, 1e-3);
    assetIsMatrixesSimilar(corrMatrix,cov1.correlationMatrix().matrix, 1e-3);

    CovarianceMatrix cov2(CorrelationMatrix{corrMatrix},variances);
    assert(cov2.correlationMatrix().matrix.rows()==3);
    assert(isClose(cov2.matrix(0,1), -0.000005150756644, 1e-9));
    assert(cov2.variance(1) == 0.00144691174);
    assetIsMatrixesSimilar(covMatrix,cov2.matrix, 1e-3);
    assetIsMatrixesSimilar(corrMatrix,cov2.correlationMatrix().matrix, 1e-3);
}

void testCovarianceMatrixNegativeVarianceError() {
    using namespace numcpp::stats; 
    Eigen::MatrixXd corrMatrix(3,3);
    corrMatrix << 
        1.0, -0.006748268411, -0.2459694651, 
        -0.006748268411, 1.0, 0.03689723423, 
        -0.2459694651, 0.03689723423, 1.0;

    Eigen::VectorXd variances(3);
    variances << 0.0004026385246,0.00144691174,-0.0005348029993;
    CovarianceMatrix cov1{CorrelationMatrix{corrMatrix}, variances};
    assert(cov1.isValid()==false);

}

void testInvalidCorrelationMatrix() {

    using namespace numcpp::stats; 
    Eigen::MatrixXd corrMatrix(3,3);
    corrMatrix << 1.0, -0.006748268411, -0.245969, 
                -0.006748268411, 1.0, 0.03689723423, 
                -0.245969, 0.98, .03689723423;
    assert(CorrelationMatrix{corrMatrix}.isValid()==false);

}

int main() {
    testCorrelationMatrix3x3(); 
    testCorrelationMatrix2x2(); 
    testCorrelationMatrix3x3MatrixConstructor(); 
    testCorrelationMatrixNonSquareError(); 
    testCorrelationMatrixNonSymmetricError(); 
    testCorrelationMatrixNonPositiveSemiDefiniteError(); 
    testCorrelationMatrixInvalidCoefficients(); 
    testCovarianceMatrixConstructor(); 
    testCovarianceMatrixNegativeVarianceError(); 
    testInvalidCorrelationMatrix();
    return 0; 
}