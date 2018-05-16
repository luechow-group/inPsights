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
#include "ExpansionSettings.h"

class RadialGaussianBasis{
public:

    explicit RadialGaussianBasis();
    
    double operator()(double r, unsigned n) const;

    double computeCoefficient(unsigned n, unsigned l, double centerToNeighborDistance, double neighborSigma) const;

    double sigmaBasisFunction(unsigned n){ return basis_[n-1].sigma(); };
    
private:
    std::vector<Gaussian> createBasis();

    Eigen::MatrixXd Smatrix() const;

    Eigen::MatrixXd radialTransform() const;

    Eigen::MatrixXd Sab(unsigned nmax) const;

    Eigen::MatrixXd calculateRadialTransform(const Eigen::MatrixXd &Sab);

    double calculateIntegral(double ai, double ri,unsigned l, double rho_ik,double beta_ik) const;
    
    std::vector<Gaussian> basis_;
    Eigen::MatrixXd Sab_, radialTransform_;
};

#endif //AMOLQCPP_RADIALGAUSSIANBASIS_H
