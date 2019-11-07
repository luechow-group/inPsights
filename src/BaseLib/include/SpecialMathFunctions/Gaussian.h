/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_RADIALGAUSSIAN_H
#define INPSIGHTS_RADIALGAUSSIAN_H

#include <Eigen/Core>
#include <utility>

class IGaussian{

public:
    double alpha() const { return alpha_; };

    double sigma() const { return sigma_; };

    double getNormalizationConstant(){ return normalizationConstant_; };

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

    double value(const Eigen::Vector3d& r) const;

    double normalizedValue(const Eigen::Vector3d& r) const;

    const Eigen::Vector3d&  center() const;


private:
    Eigen::Vector3d rCenter_;
};

#endif //INPSIGHTS_RADIALGAUSSIAN_H
