//
// Created by Michael Heuer on 11.03.18.
//

#ifndef AMOLQCPP_RADIALBASIS_H
#define AMOLQCPP_RADIALBASIS_H

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

class RadialBasis{
public:
    RadialBasis(int nmax = 4, double rCutoff = 2.0);

    double NormalizationConstant(double rCutoff, double alpha) const;

    double phi(double r,double rCutoff, double alpha) const;

    Eigen::MatrixXd Sab(int nmax) const;

    Eigen::MatrixXd inverseMatrixSqrt(const Eigen::MatrixXd& mat) const;

    double operator()(double r, int idx) const;

    int nmax() const;

private:
    double rCutoff_;
    Eigen::MatrixXd radialTransform_;
};

#endif //AMOLQCPP_RADIALBASIS_H
