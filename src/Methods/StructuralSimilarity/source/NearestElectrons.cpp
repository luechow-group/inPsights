//
// Created by leonard on 13.05.19.
//

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
                Elements::ElementInfo::Z(elementType) - Elements::ElementInfo::valElectrons(elementType);
        return getNearestElectronsIndices(electrons, nucleus.position(), coreElectronsNumber);
    };

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei) {
        // returns all core electron indices of a given vector of nuclei
        std::list<long> indices;

        for (int k = 0; k < nuclei.numberOfEntities(); ++k)
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
        std::priority_queue<std::pair<double, long>> q;
        for (long i = 0; i < electrons.numberOfEntities(); i++) {
            q.push(std::pair<double, long>(-distanceFunction(electrons[i].position(), positions), i));
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
        std::priority_queue<std::pair<double, int>> q;
        for (long i = 0; i < electrons.numberOfEntities(); i++) {
            q.push(std::pair<double, long>(-Metrics::distance(electrons[i].position(), position),
                                           i));
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
