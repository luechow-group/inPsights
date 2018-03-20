//
// Created by Michael Heuer on 20.03.18.
//

#ifndef AMOLQCPP_BOOSTSPHERICALHARMONICS_H
#define AMOLQCPP_BOOSTSPHERICALHARMONICS_H

#include <boost/math/special_functions/spherical_harmonic.hpp>

class BoostSphericalHarmonics{
public:

    static double realSphericalHarmonicY(unsigned l, int m, double theta, double phi){
        using namespace std;
        using namespace boost::math;
        if ( m < 0) {
            return sqrt(2)*pow(-1,m)* spherical_harmonic_i<double>(l,abs(m), theta, phi);
        } else if (m > 0){
            return sqrt(2)*pow(-1,m)* spherical_harmonic_r<double>(l, m, theta, phi);
        } else {
            return spherical_harmonic_r<double>(l, 0, theta, phi);
        }
    }
};

#endif //AMOLQCPP_BOOSTSPHERICALHARMONICS_H
