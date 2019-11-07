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

#include <NearestElectrons.h>
#include <Metrics.h>
#include <ElementInfo.h>
#include <queue>
#include <algorithm>

namespace NearestElectrons {
    std::list<long>
    getNonValenceIndices(const ElectronsVector &electrons, const Atom &nucleus) {
        // returns the core electron indices of a given nucleus
        const Elements::ElementType &elementType = nucleus.type();
        const long coreElectronsNumber =
                Elements::ElementInfo::Z(elementType) - Elements::ElementInfo::valenceElectrons(elementType);
        return getNearestElectronsIndices(electrons, nucleus.position(), coreElectronsNumber);
    };

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei) {
        // returns all core electron indices of a given vector of nuclei
        std::list<long> indices;

        for (long k = 0; k < nuclei.numberOfEntities(); ++k)
            indices.splice(indices.end(), getNonValenceIndices(electrons, nuclei[k]));
        indices.sort();
        return indices;
    };

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const std::vector<Eigen::Vector3d> &positions,
                               const long &maximalCount, const double &maximalDistance,
                               std::function<double(const Eigen::Vector3d &,
                                                    const std::vector<Eigen::Vector3d> &)>
                               &distanceFunction,
                               const bool &valenceOnly) {
        // returns indices of electrons closest to a vector of 'positions'
        // the metric is defined by the 'distanceFunction'
        // the indices are chosen beginning with the closest one, until 'maximalCount' or 'maximalDistance' is reached
        std::priority_queue<
                std::pair<double, int>,
                std::vector<std::pair<double, int>>,
                std::greater<std::pair<double, int>>> q;

        for (long i = 0; i < electrons.numberOfEntities(); i++) {
            q.push({distanceFunction(electrons[i].position(), positions), i});
        };

        std::list<long> excludedIndices;
        if (valenceOnly) excludedIndices = getNonValenceIndices(electrons, nuclei);

        std::list<long> indices;
        for (long j = 0; j < maximalCount; ++j) {
            // pop queue if index in pair q.top() is in excludedIndices
            while (std::find(excludedIndices.begin(), excludedIndices.end(), q.top().second) != excludedIndices.end()) {
                q.pop();
            };

            if (q.top().first > maximalDistance) break;
            indices.emplace_back(q.top().second);
            q.pop();
        };

        return indices;
    };

    std::list<long> getNearestElectronsIndices(const ElectronsVector &electrons,
                                               const Eigen::Vector3d &position,
                                               const long &count) {
        // returns indices of the 'count' electrons closest to 'position' (without restrictions)
        std::priority_queue<
        std::pair<double, int>,
        std::vector<std::pair<double, int>>,
        std::greater<std::pair<double, int>>> q;

        for (long i = 0; i < electrons.numberOfEntities(); i++) {
            q.push({Metrics::distance(electrons[i].position(), position),i});
        };

        std::list<long> indices;
        for (int i = 0; i < count; ++i) {
            indices.emplace_back(q.top().second);
            q.pop();
        };

        indices.sort();

        return indices;
    };
}