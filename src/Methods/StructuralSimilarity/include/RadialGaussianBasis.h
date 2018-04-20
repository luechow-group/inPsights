//
// Created by Michael Heuer on 22.03.18.
//

#ifndef AMOLQCPP_RADIALGAUSSIANBASIS_H
#define AMOLQCPP_RADIALGAUSSIANBASIS_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions.hpp>
#include "Gaussian.h"

class RadialGaussianBasis{
public:

    // adaptive
    explicit RadialGaussianBasis(int nmax, int lmax, double sigma0 = 1/2.);

    // equispaced
    explicit RadialGaussianBasis(int nmax, double rCut = 4.0, double sigma = 1/2.);

    double operator()(double r, int n) const;

private:
    std::vector<Gaussian> createBasis(int nmax, double rCut, double sigma);

    std::vector<Gaussian> createBasis(int nmax, int lmax, double sigma0);

    Eigen::MatrixXd Smatrix() const;

    Eigen::MatrixXd radialTransform() const;

    Eigen::MatrixXd Sab(int nmax) const;

    Eigen::MatrixXd calculateRadialTransform(const Eigen::MatrixXd &Sab);

    int nmax_;
    double sigma0_;
    std::vector<Gaussian> basis_;
    double rCut_;
    Eigen::MatrixXd Sab_, radialTransform_;
};

#endif //AMOLQCPP_RADIALGAUSSIANBASIS_H
