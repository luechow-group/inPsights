//
// Created by Michael Heuer on 20.03.18.
//

#ifndef INPSIGHTS_BOOSTSPHERICALHARMONICS_H
#define INPSIGHTS_BOOSTSPHERICALHARMONICS_H

#include <Eigen/Core>

namespace BoostSphericalHarmonics{

    Eigen::Vector3d ToDirection(double phi, double theta);

    void toSphericalCoords(const Eigen::Vector3d &dir, double &theta, double &phi);

    void toSphericalCoordsStandardizedWith2PiShift(const Eigen::Vector3d &vec, double &r, double &theta, double &phi);

    double realSphericalHarmonicY(unsigned l, int m, double theta, double phi);

    double realSphericalHarmonicY(unsigned l, int m, const Eigen::Vector3d& dir);

    // Clamp the first argument to be greater than or equal to the second
    // and less than or equal to the third.
    double Clamp(double val, double min, double max);

    // Return true if the first value is within epsilon of the second value.
    bool NearByMargin(double actual, double expected);

} // namespace BoostSphericalHarmonics

#endif //INPSIGHTS_BOOSTSPHERICALHARMONICS_H
