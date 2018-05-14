//
// Created by Michael Heuer on 20.03.18.
//

#include <BoostSphericalHarmonics.h>

using namespace std;
using namespace boost::math;

namespace BoostSphericalHarmonics {

    Eigen::Vector3d ToVector(double phi, double theta) {
        double r = sin(theta);
        return {r * cos(phi), r * sin(phi), cos(theta)};
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


    void ToSphericalCoords(const Eigen::Vector3d &dir, double &theta, double &phi) {
        assert(NearByMargin(dir.squaredNorm(), 1.0) && "dir is not unit");
        // Explicitly clamp the z coordinate so that numeric errors don't cause it
        // to fall just outside of acos' domain.
        theta = acos(Clamp(dir.z(), -1.0, 1.0));
        // We don't need to divide dir.y() or dir.x() by sin(theta) since they are
        // both scaled by it and atan2 will handle it appropriately.
        phi = atan2(dir.y(), dir.x());
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
        ToSphericalCoords(dir, theta, phi);

        return realSphericalHarmonicY(l, m, theta, phi);
    }
} // namespace BoostSphericalHarmonics