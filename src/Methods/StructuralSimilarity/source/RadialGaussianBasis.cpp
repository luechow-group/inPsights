//
// Created by Michael Heuer on 20.04.18.
//

#include "RadialGaussianBasis.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <unsupported/Eigen/MatrixFunctions>

// adaptive
RadialGaussianBasis::RadialGaussianBasis(unsigned nmax, unsigned lmax, double sigma0)
        : nmax_(nmax),
          lmax_(lmax),
          sigma0_(sigma0),
          basis_(createBasis(nmax,lmax,sigma0)),
          rCut_((*basis_.end()).center()),//TODO check this!
          Sab_(Sab(nmax)),
          radialTransform_(calculateRadialTransform(Sab_))
{};

// equispaced
RadialGaussianBasis::RadialGaussianBasis(unsigned nmax, double rCut, unsigned lmax, double sigma)
        : nmax_(nmax),
          lmax_(lmax),
          sigma0_(sigma),
          basis_(createBasis(nmax,rCut,sigma)),
          rCut_(rCut),
          Sab_(Sab(nmax)),
          radialTransform_(calculateRadialTransform(Sab_))
{};

double RadialGaussianBasis::operator()(double r, unsigned n) const{
    assert(n > 0 && "The radial basis function index must be positive");
    assert(n <= nmax_ && "The radial basis function index must be smaller than or equal to nmax");

    Eigen::VectorXd hvec(nmax_);

    for (int i = 0; i < nmax_; ++i) {
        hvec(i) = basis_[i].g2_r2_normalizedValue(r); //TODO store this?
    }
    return radialTransform_.col(n-1).dot(hvec);
};

Eigen::MatrixXd RadialGaussianBasis::Smatrix() const{
    return Sab_;
};

Eigen::MatrixXd RadialGaussianBasis::radialTransform() const{
    return radialTransform_;
};

std::vector<Gaussian> RadialGaussianBasis::createBasis(unsigned nmax, double rCut, double sigma) {

    std::vector<Gaussian> basis;
    for (int i = 0; i < nmax; ++i) {
        double rCenter = (rCut*i)/double(nmax);
        basis.emplace_back(rCenter,sigma);
    }
    return basis;
};

std::vector<Gaussian> RadialGaussianBasis::createBasis(unsigned nmax, unsigned lmax, double sigma0) {
    double rCenter = 0;
    double sigmaStride = 1/2.;

    std::vector<Gaussian> basis;
    for (int i = 1; i <= nmax; ++i) {
        double sigma = std::sqrt( 4/(2.*lmax + 1) * rCenter*rCenter + sigma0*sigma0 );
        basis.emplace_back(rCenter,sigma);
        rCenter += sigmaStride*sigma;
    }
    return basis;
};

Eigen::MatrixXd RadialGaussianBasis::Sab(unsigned nmax) const{
    Eigen::MatrixXd S(nmax,nmax);

    for (int i = 0; i < basis_.size(); ++i) {
        for (int j = 0; j < basis_.size(); ++j) { // skip iterations
            double a = basis_[i].alpha();
            double b = basis_[j].alpha();
            double rCenterA = basis_[i].center();
            double rCenterB = basis_[j].center();

            double w = a+b;
            double W0 = a*rCenterA + b*rCenterB;
            double s;
            s = 1./(4.*pow(w, 2.5));
            s *= exp(-a*rCenterA*rCenterA-b*rCenterB*rCenterB);
            s *= 2*sqrt(w)*W0
                 + sqrt(M_PI)*exp(std::pow(W0,2)/w)*(w+2*std::pow(W0,2))
                   *boost::math::erfc<double>(-W0/sqrt(w));
            s *= basis_[i].normalizationConstant_g2_r2()*basis_[j].normalizationConstant_g2_r2();
            S(i,j) = s;
        }
    }
    return S;
};

Eigen::MatrixXd RadialGaussianBasis::calculateRadialTransform(const Eigen::MatrixXd &Sab){
    Eigen::LLT<Eigen::MatrixXd> lltOfSab(Sab);
    auto L = lltOfSab.matrixL();
    return Eigen::Inverse<Eigen::MatrixXd>(L);
}
