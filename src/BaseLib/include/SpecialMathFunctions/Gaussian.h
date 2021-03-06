// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_RADIALGAUSSIAN_H
#define INPSIGHTS_RADIALGAUSSIAN_H

#include <Eigen/Core>
#include <utility>
#include <NaturalConstants.h>

class IGaussian{

public:
    double alpha() const { return alpha_; };

    double sigma() const { return sigma_; };

    double getNormalizationConstant(){ return normalizationConstant_; };

    static double calculateAlpha(double sigma);

private:
    virtual double calculateNormalizationConstant() const = 0;

protected:
    explicit IGaussian(double sigma = 1/2.);

    double sigma_,alpha_;
    double normalizationConstant_;
};


class Gaussian : public IGaussian{
public:
    explicit Gaussian(double rCenter = 0, double sigma = 1/2.);

    // computes the normalization constant for the 1D volume integral
    double calculateNormalizationConstant() const override;

    // computes the normalization constant for the integral S g^2 r^2 dr
    double normalizationConstant_g2_r2() const;

    // computes the normalization constant for the integral S g r^2 dr
    double normalizationConstant_g_r2() const;

    double value(double r) const;

    double g2_r2_normalizedValue(double r) const;

    double center() const;

private:
    double rCenter_;
    double normalizationConstant_g2_r2_;
};


class SphericalGaussian : public IGaussian{
public:
    SphericalGaussian(const Eigen::Vector3d& rCenter = Eigen::Vector3d::Zero(), double sigma = 1/2.);

    double calculateNormalizationConstant() const override;

    // computes the normalization constant for the 3D volume integral without the need of creating an object
    static double calculateNormalizationConstant(double sigma) {
        return std::pow(calculateAlpha(sigma)/Constant::pi, 3./2.);
    }

    double value(const Eigen::Vector3d& r) const;

    double normalizedValue(const Eigen::Vector3d& r) const;

    const Eigen::Vector3d&  center() const;


private:
    Eigen::Vector3d rCenter_;
};

#endif //INPSIGHTS_RADIALGAUSSIAN_H
