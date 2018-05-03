//
// Created by Michael Heuer on 18.04.18.
//

#ifndef AMOLQCPP_ANGULARBASIS_H
#define AMOLQCPP_ANGULARBASIS_H

#include <Eigen/Core>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <BoostSphericalHarmonics.h>
#include "RadialGaussianBasis.h"

class AngularBasis{
public:
    explicit AngularBasis(unsigned lmax);

    std::complex<double> computeCoefficient(unsigned l, int m, const Eigen::Vector3d& position) const;

    std::complex<double> computeCoefficient(unsigned l, int m, double theta, double phi) const;

private:
    unsigned lmax_;
};

#endif //AMOLQCPP_ANGULARBASIS_H
