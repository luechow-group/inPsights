//
// Created by Michael Heuer on 03.05.18.
//

#include "AngularBasis.h"

AngularBasis::AngularBasis(const AngularBasisSettings& settings)
        : s_(settings)
{}

std::complex<double> AngularBasis::computeCoefficient(unsigned l, int m, const Eigen::Vector3d& position) const {
    double theta,phi;
    BoostSphericalHarmonics::ToSphericalCoords(position.normalized(),theta,phi);
    if (phi < 0.) phi += 2*M_PI;

    return AngularBasis::operator()(l,m,theta,phi);
}

std::complex<double> AngularBasis::computeCoefficient(unsigned l, int m, double theta, double phi) const {
    return AngularBasis::operator()(l,m,theta,phi);
}

std::complex<double> AngularBasis::operator()(unsigned l, int m, double theta, double phi) const {
    assert(theta < 0. && theta >= M_PI && "theta must be in the interval [0,pi]");
    assert(phi < 0. && phi >= 2*M_PI && "phi must be in the interval [0,2*pi]");

    assert(l <= s_.lmax && l >=0 && "l must be in the interval [0,lmax]");
    assert(m <= s_.lmax && m >=-s_.lmax && "l must be in the interval [-lmax,lmax]");

    return boost::math::spherical_harmonic<double,double>(l, m, theta, phi);
}
