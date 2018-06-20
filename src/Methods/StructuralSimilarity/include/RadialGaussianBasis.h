//
// Created by Michael Heuer on 22.03.18.
//

#ifndef AMOLQCPP_RADIALGAUSSIANBASIS_H
#define AMOLQCPP_RADIALGAUSSIANBASIS_H

#include "Gaussian.h"
#include <Eigen/Core>
#include "ExpansionSettings.h"
#include <vector>

struct ModifiedSphericalBessel1stKind
{
    ModifiedSphericalBessel1stKind(int degree);
    void evaluate(double r, bool differentiate);

    static std::vector<double> eval(int degree, double r);
    static constexpr double RADZERO = 1e-10;
    static constexpr double SPHZERO = 1e-4;

    int _degree;
    std::vector<double> _in;
    std::vector<double> _din;
};

class RadialGaussianBasis{
public:

    explicit RadialGaussianBasis();
    
    double operator()(double r, unsigned n) const;

    Eigen::MatrixXd computeCoefficients(double centerToNeighborDistance, double neighborSigma) const;

    double sigmaBasisFunction(unsigned n){ return basis_[n-1].sigma(); };

private:
    std::vector<Gaussian> createBasis();

    Eigen::MatrixXd Smatrix() const;

    Eigen::MatrixXd radialTransform() const;

    Eigen::MatrixXd Sab(unsigned nmax) const;

    Eigen::MatrixXd calculateRadialTransform(const Eigen::MatrixXd &Sab);

    std::vector<double> calculateIntegrals(double ai, double ri, double rho_ik,double beta_ik) const;
    
    std::vector<Gaussian> basis_;
    Eigen::MatrixXd Sab_, radialTransform_;
};

#endif //AMOLQCPP_RADIALGAUSSIANBASIS_H
