//
// Created by Michael Heuer on 18.04.18.
//

#ifndef AMOLQCPP_ANGULARBASIS_H
#define AMOLQCPP_ANGULARBASIS_H

#include <Eigen/Core>
#include <boost/math/special_functions/spherical_harmonic.hpp>

static const double RADZERO = 1e-10;

class AngularBasis{
public:

    AngularBasis(unsigned lmax = 4)
            : lmax_(lmax)
    {}

    std::complex<double> operator()(unsigned l, int m, const Eigen::Vector3d& r) const {
        double theta = acos(r.z());
        double phi = atan2(r.y(), r.x());
        if (phi < 0.) phi += 2*M_PI;
        this->(l,m,theta,phi);
    };

    std::complex<double> operator()(unsigned l, int m, double theta, double phi) const {
        assert(theta < 0. && theta >= M_PI && "theta must be in the interval [0,pi]");
        assert(phi < 0. && phi >= 2*M_PI && "phi must be in the interval [0,2*pi]");

        assert(l <= lmax_ && l >=0 && "l must be in the interval [0,lmax]");
        assert(m <= lmax_ && m >=-lmax_ && "l must be in the interval [-lmax,lmax]");

        return boost::math::spherical_harmonic<double,double>(l, m, theta, phi);
    };

private:
    unsigned lmax_;
    Eigen::VectorXcd coeffs;
};

#endif //AMOLQCPP_ANGULARBASIS_H
