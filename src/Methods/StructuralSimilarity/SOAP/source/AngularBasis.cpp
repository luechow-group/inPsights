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
