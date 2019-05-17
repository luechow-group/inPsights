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
    getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &k) {
        const Elements::ElementType &elementType = nuclei[k].type();
        const long coreElectronsNumber =
                Elements::ElementInfo::Z(elementType) - Elements::ElementInfo::valElectrons(elementType);
        return getNearestElectronsIndices(electrons, nuclei, k, coreElectronsNumber);
    };

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei) {
        std::list<long> indices;

        for (int k = 0; k < nuclei.numberOfEntities(); ++k)
            indices.splice(indices.end(), getNonValenceIndices(electrons, nuclei, k));
        indices.sort();
        return indices;
    };

    ElectronsVector
    getNearestValenceElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const Eigen::Vector3d &position,
                               const long &count) {
        return electrons[getNearestValenceIndices(electrons, nuclei, position, count)];
    };

    ElectronsVector
    getNearestValenceElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index,
                               const long &count) {
        return getNearestValenceElectrons(electrons, nuclei, nuclei[index].position(), count);
    };

    ElectronsVector
    getNearestValenceElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index1,
                               const int &index2,
                               const long &count) {
        Eigen::Vector3d position =
                (nuclei[index1].position() + nuclei[index2].position()) / 2;
        return getNearestValenceElectrons(electrons, nuclei, position, count);
    };

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                             const Eigen::Vector3d &position,
                             const long &count) {
        std::priority_queue<std::pair<double, long>> q;
        for (long i = 0; i < electrons.numberOfEntities(); i++) {
            q.push(std::pair<double, long>(-Metrics::distance(electrons[i].position(), position),
                                           i));
        };

        std::list<long> excludedIndices = getNonValenceIndices(electrons, nuclei);

        std::list<long> indices;
        for (long j = 0; j < count; ++j) {
            // pop queue if index in pair q.top() is in excludedIndices
            while (std::find(excludedIndices.begin(), excludedIndices.end(), q.top().second) != excludedIndices.end()) {
                q.pop();
            };
            indices.emplace_back(q.top().second);
            q.pop();
        };

        indices.sort();

        return indices;
    };

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index,
                             const long &count) {
        return getNearestValenceIndices(electrons, nuclei, nuclei[index].position(),
                                        count);
    };

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index1,
                             const int &index2,
                             const long &count) {
        Eigen::Vector3d position =
                (nuclei[index1].position() + nuclei[index2].position()) / 2;
        return getNearestValenceIndices(electrons, nuclei, position, count);
    };

    std::list<long> getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                                               const Eigen::Vector3d &position,
                                               const long &count) {
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

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index,
                               const long &count) {
        return getNearestElectronsIndices(electrons, nuclei, nuclei[index].position(),
                                          count);
    };

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index1,
                               const int &index2,
                               const long &count) {
        Eigen::Vector3d position =
                (nuclei[index1].position() + nuclei[index2].position()) / 2;
        return getNearestElectronsIndices(electrons, nuclei, position, count);
    };
}
