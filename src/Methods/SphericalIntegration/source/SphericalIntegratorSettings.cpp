//
// Created by Michael Heuer on 14.03.18.
//

#include "SphericalIntegratorSettings.h"

SphericalIntegratorSettings SphericalIntegratorSettings::defaults() {
    SphericalIntegratorSettings s{};
    s.lebedevRule = Lebedev::OrderType::LD0006;

    s.radialGridPoints = 100;
    s.gaussKronrodRule = Eigen::Integrator<double>::QuadratureRule::GaussKronrod15;
    s.desiredAbsoluteRadialIntegrationError = 0.0;
    s.desiredRelativeRadialIntegrationError = Eigen::NumTraits<double>::epsilon()*50;
    return s;
}

/*
 * Automatic settings for expansion in spherical harmonics (with maximal polynomial degree lmax)
 * and radial basis (with maximal polynomial degree nmax) for a given cutoff radius rCutoff.
 */
SphericalIntegratorSettings SphericalIntegratorSettings::expansion(unsigned nmax, unsigned lmax, double rCutoff) {
    SphericalIntegratorSettings s{};
    s.lebedevRule = Lebedev::findAdequateRule(lmax);
    assert(s.lebedevRule != Lebedev::OrderType::NotAvailable && "Lmax can not be too large.");
    assert(rCutoff > 0 && "rCutoff must be positive.");

    // The radial basis we employ is a polynomial expansion with order nmax+2
    s.gaussKronrodRule = GaussKronrod::findAdequateRule(nmax + 2);
    s.radialGridPoints = unsigned(std::round(100*std::abs(rCutoff)));

    s.desiredAbsoluteRadialIntegrationError = 0.0;
    s.desiredRelativeRadialIntegrationError = Eigen::NumTraits<double>::epsilon()*50;

    return s;
}