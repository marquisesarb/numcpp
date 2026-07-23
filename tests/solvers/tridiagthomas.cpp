#include <cassert>
#include <numcpp/solvers/tridiagthomas.hpp>

bool isClose(double a, double b, double tol = 1e-3) {return std::abs(a - b) <= tol;}

void testThomas() {

    //exemple here: https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://matlabgeeks.weebly.com/uploads/8/0/4/8/8048228/thomas_algorithm_and_tridiagonal_matrix-_example.pdf&ved=2ahUKEwix68fOh8mVAxW41QIHHXYqMd8QFnoECCgQAQ&usg=AOvVaw1pJ3WcPX-sjJCsva3hXC-C

    Eigen::SparseMatrix<double> M(4, 4);
    std::vector<Eigen::Triplet<double>> coefficients;
    coefficients.push_back({Eigen::Triplet<double>(0,0,-2.6)});
    coefficients.push_back({Eigen::Triplet<double>(0,1,1.0)});
    coefficients.push_back({Eigen::Triplet<double>(3,2,1.0)});
    coefficients.push_back({Eigen::Triplet<double>(3,3,-2.6)});
    for (size_t i = 1; i<3; i++) {

        coefficients.push_back({Eigen::Triplet<double>(i,i-1,1.0)});
        coefficients.push_back({Eigen::Triplet<double>(i,i,-2.6)});
        coefficients.push_back({Eigen::Triplet<double>(i,i+1,1.0)});
    }
    M.setFromTriplets(coefficients.begin(), coefficients.end());
    Eigen::VectorXd d(4); 
    d << -240.0, 0.0, 0.0, -150.0; 

    Eigen::VectorXd d2 = M*numcpp::solvers::tridiagonalMatrixSolverThomas(M, d);

    assert(isClose(d(0), d2(0), 1e-10));
    assert(isClose(d(1), d2(1), 1e-10));
    assert(isClose(d(2), d2(2), 1e-10));
    assert(isClose(d(3), d2(3), 1e-10));

}

void testThomasCyclic() {

    Eigen::SparseMatrix<double> M(4, 4);
    std::vector<Eigen::Triplet<double>> coefficients;
    coefficients.push_back({Eigen::Triplet<double>(0,3,1.0)});
    coefficients.push_back({Eigen::Triplet<double>(3,0,1.0)});
    coefficients.push_back({Eigen::Triplet<double>(0,0,4.0)});
    coefficients.push_back({Eigen::Triplet<double>(0,1,1.0)});
    coefficients.push_back({Eigen::Triplet<double>(3,2,1.0)});
    coefficients.push_back({Eigen::Triplet<double>(3,3,4.0)});
    for (size_t i = 1; i<3; i++) {

        coefficients.push_back({Eigen::Triplet<double>(i,i-1,1.0)});
        coefficients.push_back({Eigen::Triplet<double>(i,i,4.0)});
        coefficients.push_back({Eigen::Triplet<double>(i,i+1,1.0)});
    }
    M.setFromTriplets(coefficients.begin(), coefficients.end());
    Eigen::VectorXd d(4); 
    d << 10.0, 12.0, 18.0, 20.0; 
    Eigen::VectorXd u(4);
    u << 1.0, 0.0, 0.0, 1.0; 

    Eigen::VectorXd d2 = M*numcpp::solvers::tridiagonalMatrixSolverThomasShermanMorrison(M, d,u,u);

    assert(isClose(d(0), d2(0), 1e-10));
    assert(isClose(d(1), d2(1), 1e-10));
    assert(isClose(d(2), d2(2), 1e-10));
    assert(isClose(d(3), d2(3), 1e-10));

}


int main() {

    testThomas();
    testThomasCyclic(); 

}