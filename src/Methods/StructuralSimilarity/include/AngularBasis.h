//
// Created by Michael Heuer on 18.04.18.
//

#ifndef AMOLQCPP_ANGULARBASIS_H
#define AMOLQCPP_ANGULARBASIS_H

#include <Eigen/Core>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <BoostSphericalHarmonics.h>
#include "RadialGaussianBasis.h"

class AngularBasis{
public:

    AngularBasis(unsigned lmax = 4)
            : lmax_(lmax),
              coefficientsYlm_(lmax_*lmax_+2*lmax_+1)
    {}

    std::complex<double> operator()(unsigned l, int m, const Eigen::Vector3d& r) const {
        double theta,phi;
        BoostSphericalHarmonics::ToSphericalCoords(r.normalized(),theta,phi);
        if (phi < 0.) phi += 2*M_PI;

        return AngularBasis::operator()(l,m,theta,phi);
    };

    std::complex<double> operator()(unsigned l, int m, double theta, double phi) const {
        assert(theta < 0. && theta >= M_PI && "theta must be in the interval [0,pi]");
        assert(phi < 0. && phi >= 2*M_PI && "phi must be in the interval [0,2*pi]");

        assert(l <= lmax_ && l >=0 && "l must be in the interval [0,lmax]");
        assert(m <= lmax_ && m >=-lmax_ && "l must be in the interval [-lmax,lmax]");

        return boost::math::spherical_harmonic<double,double>(l, m, theta, phi);
    };


    void computeCoefficients(const Eigen::Vector3d& rvec) {

        double r = rvec.norm();

        if (r < ZeroLimits::radiusZero) {
            coefficientsYlm_(0) = boost::math::spherical_harmonic<double, double>(0, 0, 0.0, 0.0);
        } else {
            double theta,phi;
            BoostSphericalHarmonics::ToSphericalCoords(rvec.normalized(),theta,phi);
            if (phi < 0.) phi += 2*M_PI;

            for (unsigned l = 0; l <= lmax_; ++l) {
                for (int m = -l; m <= l; ++m) {
                    coefficientsYlm_(l*l+l+m) = boost::math::spherical_harmonic<double, double>(l,m,theta,phi);
                }
            }
        }
    }

    std::complex<double> getCoefficient(unsigned l, int m){
        return coefficientsYlm_(l*l+l+m);
    }

private:
    unsigned lmax_;
    Eigen::VectorXcd coefficientsYlm_; //TODO why complex?
};

#endif //AMOLQCPP_ANGULARBASIS_H
