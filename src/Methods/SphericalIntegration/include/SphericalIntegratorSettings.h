//
// Created by Michael Heuer on 14.03.18.
//

#ifndef AMOLQCPP_SPHERICALINTEGRATORSETTINGS_H
#define AMOLQCPP_SPHERICALINTEGRATORSETTINGS_H

#include "GaussKronrodCartesianIntegration.h"
#include "LebedevSphericalIntegration/OrderType.h"

class SphericalIntegratorSettings {
public:
    Lebedev::OrderType lebedevRule;
    unsigned radialGridPoints;
    GaussKronrod::Integrator<double>::QuadratureRule gaussKronrodRule;
    double desiredRelativeRadialIntegrationError;
    double desiredAbsoluteRadialIntegrationError;

    static SphericalIntegratorSettings defaults();

    static SphericalIntegratorSettings expansion(unsigned lmax, unsigned nmax,  double rCutoff);
};

#endif //AMOLQCPP_SPHERICALINTEGRATORSETTINGS_H
