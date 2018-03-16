//
// Created by Michael Heuer on 11.03.18.
//

#include "RadialBasis.h"

RadialBasis::RadialBasis(int nmax, double rCutoff)
        : rCutoff_(rCutoff),
          W_(W(nmax)){
};

double RadialBasis::NormalizationConstant(double rCutoff, double alpha) const {
    return std::sqrt( std::pow(rCutoff, 2*alpha+5) / (2*alpha+5) );
}

double RadialBasis::phi(double r,double rCutoff, double alpha) const {
    return  std::pow( rCutoff-r, alpha+2) / NormalizationConstant(rCutoff,alpha);
}

Eigen::MatrixXd RadialBasis::Sab(int nmax) const {
    assert(nmax > 0 && "The number of radial basis functions must be greater than zero.");
    Eigen::MatrixXd Sab (nmax,nmax);

    for (int i = 1; i <= Sab.rows(); ++i) {
        for (int j = 1; j <= nmax; ++j) {
            Sab(i-1,j-1) = std::sqrt((5+2*i)*(5+2*j))/double(5+i+j);
        }
    }
    return Sab;
}

Eigen::MatrixXd RadialBasis::W(int nmax) const {
    assert(nmax > 0 && "The number of radial basis functions must be positive.");
    return Sab(nmax).inverse().sqrt();
}

double RadialBasis::operator()(double r, int idx) const {
    int nMax = nmax();
    assert(idx >= 0 && "The radial basis function index must be smaller than nmax");
    assert(idx < nMax && "The radial basis function index must be smaller than nmax");

    double sum = 0;
    for (int a = 0; a < nMax; ++a) {
        sum += W_(idx,a) * phi(r, rCutoff_,a);
    }
    return sum;
};

int RadialBasis::nmax() const {
    assert(W_.rows() == W_.cols() && "The Matrix W must be square.");
    assert(W_.rows() > 0 && "Nmax must be positive.");
    assert(W_.rows() < std::numeric_limits<int>::max() && "The dimensions of the matrix should fit into an integer.");
    return static_cast<int>(W_.rows());
};
