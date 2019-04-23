//
// Created by Michael Heuer on 03.05.18.
//

#include <AngularBasis.h>
#include <SOAPSettings.h>
#include <SpecialMathFunctions/BoostSphericalHarmonics.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace SOAP;

std::complex<double> AngularBasis::computeCoefficient(unsigned l, int m, const Eigen::Vector3d& position) {
    double theta,phi;
    BoostSphericalHarmonics::toSphericalCoords(position.normalized(), theta, phi);
    if (phi < 0.) phi += 2*M_PI;

    return AngularBasis::computeCoefficient(l,m,theta,phi);
}

std::complex<double> AngularBasis::computeCoefficient(unsigned l, int m, double theta, double phi) {
    assert(theta >= 0. && theta <= M_PI && "theta must be in the interval [0,pi]");
    assert(phi >= 0. && phi <= 2*M_PI && "phi must be in the interval [0,2*pi]");

    assert(l <= Angular::settings.lmax() && l >=0 && "l must be in the interval [0,lmax]");
    assert(unsigned(abs(m)) <= Angular::settings.lmax() && "l must be in the interval [-lmax,lmax]");


    //TODO? use real spherical harmonic?
    //BoostSphericalHarmonics::realSphericalHarmonicY(l,m,theta,phi);
    return boost::math::spherical_harmonic<double,double>(l, m, theta, phi);
}
