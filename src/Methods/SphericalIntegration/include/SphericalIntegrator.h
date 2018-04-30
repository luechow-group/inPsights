//
// Created by Michael Heuer on 12.03.18.
//

#ifndef AMOLQCPP_SPHERICALINTEGRATOR_H
#define AMOLQCPP_SPHERICALINTEGRATOR_H

#include "SphericalIntegratorSettings.h"
#include "LebedevSphericalIntegration/SphericalSurfaceIntegrator.h"
#include "GaussKronrodCartesianIntegration.h"

class SpatialFunction;

class SphericalIntegrator {
    // Eigen::Integrator is a friend to allow for access to the private method double operator()(double r) const;
    template <typename Scalar> friend class Eigen::Integrator;

public:
    explicit SphericalIntegrator(const SphericalIntegratorSettings& settings);

    double integrate(const SpatialFunction& f, double a, double b);

    double integrate(const SpatialFunction& f, double b);

private:
    void setFunction(const SpatialFunction& f);

    double operator()(double r) const;

    SpatialFunction * fPtr_;
    SphericalIntegratorSettings settings_;
    Lebedev::SphericalSurfaceIntegrator sphericalSurfaceIntegrator_;
    Eigen::Integrator<double> cartesianIntegrator_; //TODO rename to radial integrator
    Eigen::Integrator<double>::QuadratureRule quadratureRule_;
    double desiredRelativeRadialIntegrationError_;
    double desiredAbsoluteRadialIntegrationError_;
};

#endif //AMOLQCPP_SPHERICALINTEGRATOR_H
