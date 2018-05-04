//
// Created by Michael Heuer on 18.04.18.
//

#ifndef AMOLQCPP_ANGULARBASIS_H
#define AMOLQCPP_ANGULARBASIS_H

#include "ExpansionSettings.h"
#include <BoostSphericalHarmonics.h>
#include <Eigen/Core>
#include <boost/math/special_functions/spherical_harmonic.hpp>



class AngularBasis{
public:
    explicit AngularBasis(const ExpansionSettings::AngularBasisSettings& settings = ExpansionSettings::AngularBasisSettings::defaults());

    std::complex<double> computeCoefficient(unsigned l, int m, const Eigen::Vector3d& position) const;

    std::complex<double> computeCoefficient(unsigned l, int m, double theta, double phi) const;

    std::complex<double> operator()(unsigned l, int m, double theta, double phi) const;

private:
    ExpansionSettings::AngularBasisSettings s_;
};

#endif //AMOLQCPP_ANGULARBASIS_H
