//
// Created by Michael Heuer on 18.04.18.
//

#ifndef AMOLQCPP_ANGULARBASIS_H
#define AMOLQCPP_ANGULARBASIS_H

#include "ExpansionSettings.h"
#include <BoostSphericalHarmonics.h>
#include <Eigen/Core>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace AngularBasis{
    std::complex<double> computeCoefficient(unsigned l, int m, const Eigen::Vector3d& position);

    std::complex<double> computeCoefficient(unsigned l, int m, double theta, double phi);
};

#endif //AMOLQCPP_ANGULARBASIS_H
