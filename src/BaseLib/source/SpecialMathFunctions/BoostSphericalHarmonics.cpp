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

#include "SpecialMathFunctions/BoostSphericalHarmonics.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace std;
using namespace boost::math;

namespace BoostSphericalHarmonics {

    Eigen::Vector3d ToDirection(double phi, double theta) {
        double r = sin(theta);
        return {r * cos(phi), r * sin(phi), cos(theta)};
    }

    Eigen::Vector3d ToVector(double r,double phi, double theta) {
        return r*ToDirection(phi, theta);
    }

    // Clamp the first argument to be greater than or equal to the second
    // and less than or equal to the third.
    double Clamp(double val, double min, double max) {
        if (val < min) {
            val = min;
        }
        if (val > max) {
            val = max;
        }
        return val;
    }

    // Return true if the first value is within epsilon of the second value.
    bool NearByMargin(double actual, double expected) {
        double diff = actual - expected;
        if (diff < 0.0) {
            diff = -diff;
        }
        // 5 bits of error in mantissa (source of '32 *')
        return diff < 32 * std::numeric_limits<double>::epsilon();
    }



    void toSphericalCoords(const Eigen::Vector3d &dir, double &theta, double &phi) {
        assert(NearByMargin(dir.squaredNorm(), 1.0) && "dir is not unit");
        // Explicitly clamp the z coordinate so that numeric errors don't cause it
        // to fall just outside of acos' domain.
        theta = acos(Clamp(dir.z(), -1.0, 1.0));
        // We don't need to divide dir.y() or dir.x() by sin(theta) since they are
        // both scaled by it and atan2 will handle it appropriately.
        phi = atan2(dir.y(), dir.x());

        //theta is element of [0,pi]
        //phi is element of [-pi,pi]
    }

    void toSphericalCoords(const Eigen::Vector3d &vec, double &r, double &theta, double &phi){
        r = vec.norm();
        toSphericalCoords(vec.normalized(),theta,phi);
    }


    void toSphericalCoordsStandardizedWith2PiShift(const Eigen::Vector3d &vec, double &r, double &theta, double &phi){
        r = vec.norm();
        
        if (r <= 0) {
            theta = 0.;
            phi = 0.;
        } else {
            toSphericalCoords(vec.normalized(), theta, phi); //phi in [-pi,pi]
            if( theta == 0. ){
                phi = 0;
            } else if (phi < 0.) {
                phi += 2 * M_PI;  //phi now in [0,2*pi]
            }
        }
    }

    double realSphericalHarmonicY(unsigned l, int m, double theta, double phi) {
        using namespace std;
        using namespace boost::math;

        if (m == 0)  {
            return spherical_harmonic_r<double>(l, 0, theta, phi);
        } else if (m < 0) {
            return sqrt(2) * pow(-1, m) * spherical_harmonic_i<double>(l, abs(m), theta, phi);
        } else { // (m > 0)
            return sqrt(2) * pow(-1, m) * spherical_harmonic_r<double>(l, m, theta, phi);
        } 
    }

    double realSphericalHarmonicY(unsigned l, int m, const Eigen::Vector3d &dir) {
        double theta, phi;
        toSphericalCoords(dir, theta, phi);

        return realSphericalHarmonicY(l, m, theta, phi);
    }
} // namespace BoostSphericalHarmonics