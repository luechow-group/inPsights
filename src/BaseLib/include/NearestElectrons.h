/* Copyright (C) 2019 Leonard Reuter.
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

#ifndef INPSIGHTS_NEARESTELECTRONS_H
#define INPSIGHTS_NEARESTELECTRONS_H

#include <ParticlesVector.h>
#include <functional>

namespace NearestElectrons {
    std::list<long>
    getNonValenceIndices(const ElectronsVector &electrons, const Atom &nucleus);

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei);

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const std::vector<Eigen::Vector3d> &positions,
                               const long &maximalCount, const double &maximalDistance,
                               std::function<double(const Eigen::Vector3d &,
                                                    const std::vector<Eigen::Vector3d> &)>
                               &distanceFunction,
                               const bool &valenceOnly);

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const Eigen::Vector3d &position,
                               const long &count);
    std::list<long> invertedIndices(const std::list<long>& indices, std::size_t size);
}

#endif //INPSIGHTS_NEARESTELECTRONS_H
