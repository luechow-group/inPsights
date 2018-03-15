//
// Created by Michael Heuer on 12.03.18.
//

#include "SphericalIntegrator.h"

SphericalIntegrator::SphericalIntegrator(const SphericalIntegratorSettings& settings)
        : function_(nullptr),
          settings_(settings),
          sphericalSurfaceIntegrator_(settings.lebedevRule),
          cartesianIntegrator_(settings.radialGridPoints),
          quadratureRule_(settings.gaussKronrodRule),
          desiredAbsoluteRadialIntegrationError_(settings.desiredAbsoluteRadialIntegrationError),
          desiredRelativeRadialIntegrationError_(settings.desiredRelativeRadialIntegrationError){};

void SphericalIntegrator::setFunction(const SpatialFunction& f){ // const SpatialFunction&
    function_ = const_cast<SpatialFunction*>(&f);
};

double SphericalIntegrator::integrate(const SpatialFunction& f, double a, double b){
    setFunction(f);
    assert(function_ && "The pointer to the function cannot be the nullptr.");
    return cartesianIntegrator_.quadratureAdaptive(*this,a,b,
                                                   desiredAbsoluteRadialIntegrationError_,
                                                   desiredRelativeRadialIntegrationError_,
                                                   quadratureRule_);
};

double SphericalIntegrator::integrate(const SpatialFunction& f, double b){
    integrate(f,0,b);
};

double SphericalIntegrator::operator()(double r) const {
    return sphericalSurfaceIntegrator_.integrate(function_,r) * std::pow(r,2);// spherical volume element
};
