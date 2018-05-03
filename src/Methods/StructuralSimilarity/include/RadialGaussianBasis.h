//
// Created by Michael Heuer on 22.03.18.
//

#ifndef AMOLQCPP_RADIALGAUSSIANBASIS_H
#define AMOLQCPP_RADIALGAUSSIANBASIS_H

#include "Gaussian.h"
#include "ExpansionSettings.h"
#include <Eigen/Core>
#include <boost/math/special_functions/bessel.hpp>
#include <GaussKronrodCartesianIntegration.h>
#include <cmath>

namespace ZeroLimits{
    const double radiusZero = 1e-10;
}

class RadialGaussianBasis{
public:

    explicit RadialGaussianBasis(const ExpansionSettings& settings = ExpansionSettings::defaults());
    
    double operator()(double r, unsigned n) const;

    double computeCoefficient(unsigned n, unsigned l, const Eigen::Vector3d& neighborPosition, double neighborSigma) const;

    double sigmaBasisFunction(unsigned n){ return basis_[n-1].sigma(); };
    
private:
    std::vector<Gaussian> createBasis(ExpansionSettings& settings);

    Eigen::MatrixXd Smatrix() const;

    Eigen::MatrixXd radialTransform() const;

    Eigen::MatrixXd Sab(unsigned nmax) const;

    Eigen::MatrixXd calculateRadialTransform(const Eigen::MatrixXd &Sab);

    double calculateIntegral(double ai, double ri,unsigned l, double rho_ik,double beta_ik) const;

    ExpansionSettings s_;
    std::vector<Gaussian> basis_;
    Eigen::MatrixXd Sab_, radialTransform_;
};

#endif //AMOLQCPP_RADIALGAUSSIANBASIS_H
