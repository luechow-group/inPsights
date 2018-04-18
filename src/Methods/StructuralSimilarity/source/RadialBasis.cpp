//
// Created by Michael Heuer on 11.03.18.
//

#include "RadialBasis.h"
#include <Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>
//#include <unsupported>


RadialBasis::RadialBasis(int nmax, double rCutoff)
        : rCutoff_(rCutoff),
          radialTransform_(inverseMatrixSqrt(Sab(nmax))){
};

double RadialBasis::NormalizationConstant(double rCutoff, double alpha) const {
    return std::sqrt( std::pow(rCutoff, 2.0*alpha+5.0) / (2.0*alpha+5.0) );
}

double RadialBasis::phi(double r,double rCutoff, double alpha) const {
    assert(alpha > 0 && alpha <= nmax() && "The basis function number ranges from 1 to nmax.");
    return  std::pow( rCutoff-r, alpha+2.0) / NormalizationConstant(rCutoff,alpha);
}

Eigen::MatrixXd RadialBasis::Sab(int nmax) const {
    assert(nmax > 0 && "The number of radial basis functions must be greater than zero.");
    Eigen::MatrixXd Sab (nmax,nmax);
    for (int i = 1; i <= Sab.rows(); ++i) {
        for (int j = 1; j <= nmax; ++j) {
            Sab(i-1,j-1) = std::sqrt((5.0+2.0*i)*(5.0+2.0*j))/(5.0+i+j);
        }
    }
    return Sab;
}

Eigen::MatrixXd RadialBasis::inverseMatrixSqrt(const Eigen::MatrixXd& mat) const {

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(mat, Eigen::ComputeEigenvectors);
    if (eigenSolver.info() != Eigen::Success) abort();

    //TODO check if necessary, if not delete
    // With cwiseAbs(), negative eigenvalues are eliminated
    //Eigen::VectorXd inverseSqrtEigenvalues= eigenSolver.eigenvalues().cwiseAbs().cwiseInverse().cwiseSqrt();

    Eigen::VectorXd inverseSqrtEigenvalues= eigenSolver.eigenvalues().cwiseInverse().cwiseSqrt();

    //std::cout << mat << std::endl << std::endl;
    //std::cout << (eigenSolver.eigenvectors()
    //              * inverseSqrtEigenvalues.asDiagonal())
    //             * eigenSolver.eigenvectors().transpose().inverse() << std::endl;

    return (eigenSolver.eigenvectors()
           * inverseSqrtEigenvalues.asDiagonal())
           * eigenSolver.eigenvectors().transpose();
}

double RadialBasis::operator()(double r, int idx) const {
    int nMax = nmax();
    assert(idx >= 0 && "The radial basis function index must be positive.");
    assert(idx < nMax && "The radial basis function index must be smaller than nmax");

    Eigen::VectorXd hvec(nMax);

    for (int n = 1; n <= nMax; ++n) {
        hvec(n-1) = phi(r, rCutoff_, n);
    }
    return radialTransform_.col(idx).dot(hvec);
};

int RadialBasis::nmax() const {
    assert(radialTransform_.rows() == radialTransform_.cols() && "The Matrix W must be square.");
    assert(radialTransform_.rows() > 0 && "Nmax must be positive.");
    assert(radialTransform_.rows() < std::numeric_limits<int>::max() && "The dimensions of the matrix should fit into an integer.");
    return static_cast<int>(radialTransform_.rows());
};
