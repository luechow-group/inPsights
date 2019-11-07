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

#include "SpecialMathFunctions/Gaussian.h"
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <NaturalConstants.h>

IGaussian::IGaussian(double sigma)
        : sigma_(sigma),
          alpha_(1./(2*pow(sigma,2))),
          normalizationConstant_(0) {}


Gaussian::Gaussian(double rCenter, double sigma)
        : IGaussian(sigma),
          rCenter_(rCenter),
          normalizationConstant_g2_r2_(normalizationConstant_g2_r2())
{
    normalizationConstant_= Gaussian::calculateNormalizationConstant();
}

// computes the normalization constant for the 1D volume integral
double Gaussian::calculateNormalizationConstant() const {
    return  sqrt(alpha_/Constant::pi);
};

double Gaussian::center() const { return rCenter_; };

// computes the normalization constant for the integral S g^2 r^2 dr
double Gaussian::normalizationConstant_g2_r2() const{
    double w = 2*alpha_;
    double W0 = 2*alpha_*rCenter_;
    double integral_r2_g2_dr =
            1./(4.*pow(w, 5/2.))
            *exp(-w*pow(rCenter_,2))
            *(2*sqrt(w)*W0
              + sqrt(Constant::pi )
                *exp(pow(W0,2)/w)
                *(w+2*pow(W0,2))
                *boost::math::erfc<double>(-W0/sqrt(w))
            );
    assert(integral_r2_g2_dr == integral_r2_g2_dr && "coeff cannot be NaN!");
    assert(integral_r2_g2_dr > 0);
    return  1./sqrt(integral_r2_g2_dr);
};

// computes the normalization constant for the integral S g r^2 dr
double Gaussian::normalizationConstant_g_r2() const{
    double w = alpha_;
    double W0 = alpha_*rCenter_;
    double integral_r2_g_dr =
            1./(4.*pow(w, 2.5))*exp(-w*pow(rCenter_,2))*(
                    2*sqrt(w)*W0 +
                    sqrt(Constant::pi )*exp(pow(W0,2)/w)*(w+2*pow(W0,2))*(
                            1 - boost::math::erf<double>(-W0/sqrt(w))
                    )
            );
    return  1./sqrt(integral_r2_g_dr);
};

double Gaussian::value(double r) const {
    return std::exp(-alpha_*pow(r-rCenter_,2));
};

double Gaussian::g2_r2_normalizedValue(double r) const {
    return value(r)*normalizationConstant_g2_r2_;
};


SphericalGaussian::SphericalGaussian(const Eigen::Vector3d& rCenter, double sigma)
        : IGaussian(sigma),
          rCenter_(rCenter)
{
    normalizationConstant_= SphericalGaussian::calculateNormalizationConstant();
};

// computes the normalization constant for the 3D volume integral
double SphericalGaussian::calculateNormalizationConstant() const {
    return pow(alpha_/Constant::pi , 3./2.);
};

double SphericalGaussian::value(const Eigen::Vector3d &r) const {
    return std::exp(-alpha_*(r-rCenter_).squaredNorm());
}

double SphericalGaussian::normalizedValue(const Eigen::Vector3d &r) const {
    return value(r) * normalizationConstant_;
}

const Eigen::Vector3d &SphericalGaussian::center() const { return rCenter_; };
