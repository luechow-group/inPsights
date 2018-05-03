//
// Created by Michael Heuer on 22.03.18.
//

#ifndef AMOLQCPP_RADIALGAUSSIANBASIS_H
#define AMOLQCPP_RADIALGAUSSIANBASIS_H

#include "Gaussian.h"
#include <Eigen/Core>
#include <boost/math/special_functions/bessel.hpp>
#include <GaussKronrodCartesianIntegration.h>
#include <cmath>

namespace ZeroLimits{
    const double radiusZero = 1e-10;
}

class RadialGaussianBasis{
public:

    //TODO understand and discriminate the sigmas
    // sigma0 = 1/2. default

    // adaptive
    explicit RadialGaussianBasis(unsigned nmax, unsigned lmax, double sigma0 = 1/2.);

    // equispaced
    explicit RadialGaussianBasis(unsigned nmax, double rCut = 4.0, unsigned lmax = 4, double sigma = 1/2.); // sensible default?

    double operator()(double r, unsigned n) const;

    double sigma(){ return sigma0_; };

    double computeCoefficient(unsigned n, unsigned l, const Eigen::Vector3d& neighborPosition, double neighborSigma) const;

private:
    std::vector<Gaussian> createBasis(unsigned nmax, double rCut, double sigma);

    std::vector<Gaussian> createBasis(unsigned nmax, unsigned lmax, double sigma0);

    Eigen::MatrixXd Smatrix() const;

    Eigen::MatrixXd radialTransform() const;

    Eigen::MatrixXd Sab(unsigned nmax) const;

    Eigen::MatrixXd calculateRadialTransform(const Eigen::MatrixXd &Sab);

    double calculateIntegral(double ai, double ri,unsigned l, double rho_ik,double beta_ik) const;

    unsigned nmax_,lmax_;
    double sigma0_;
    std::vector<Gaussian> basis_;
    double rCut_;

    Eigen::MatrixXd Sab_, radialTransform_;
};

#endif //AMOLQCPP_RADIALGAUSSIANBASIS_H
