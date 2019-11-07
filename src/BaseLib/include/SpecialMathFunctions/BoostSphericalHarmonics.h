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
