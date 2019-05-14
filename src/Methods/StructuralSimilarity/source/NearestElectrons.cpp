//
// Created by leonard on 13.05.19.
//

#include <NearestElectrons.h>
#include <Metrics.h>
#include <ElementInfo.h>
#include <queue>
#include <algorithm>

namespace NearestElectrons {
    Eigen::ArrayXi longIndexListToPositionArrayXi(std::list<long> list) {
        Eigen::ArrayXi array(3 * list.size());
        int i = 0;
        for (auto element : list) {
            for (int j = 0; j < 3; j++) {
                array[i] = 3 * element + j;
                i++;
            };
        };
        return array;
    };

    std::list<long>
    getNonValenceIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei, const int &k) {
        const Elements::ElementType &elementType = nuclei[k].type();
        const long coreElectronsNumber =
                Elements::ElementInfo::Z(elementType) - Elements::ElementInfo::valElectrons(elementType);
        return getNearestElectronsIndices(electronPositions, nuclei, k, coreElectronsNumber);
    };

    std::list<long> getNonValenceIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei) {
        std::list<long> indices;

        for (int k = 0; k < nuclei.numberOfEntities(); ++k)
            indices.splice(indices.end(), getNonValenceIndices(electronPositions, nuclei, k));
        indices.sort();
        return indices;
    };

    Eigen::VectorXd
    getNearestValencePositions(const PositionsVector &electronPositions, const AtomsVector &nuclei,
                               const Eigen::Vector3d &position,
                               const long &count) {
        return longIndexListToPositionArrayXi(
                getNearestValenceIndices(electronPositions, nuclei, position, count)).unaryExpr(
                electronPositions.asEigenVector());
    };

    Eigen::VectorXd
    getNearestValencePositions(const PositionsVector &electronPositions, const AtomsVector &nuclei, const int &index,
                               const long &count) {
        return getNearestValencePositions(electronPositions, nuclei, nuclei[index].position(),
                                          count);
    };

    Eigen::VectorXd
    getNearestValencePositions(const PositionsVector &electronPositions, const AtomsVector &nuclei, const int &index1,
                               const int &index2,
                               const long &count) {
        Eigen::Vector3d position =
                (nuclei[index1].position() + nuclei[index2].position()) / 2;
        return getNearestValencePositions(electronPositions, nuclei, position, count);
    };

    std::list<long>
    getNearestValenceIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei,
                             const Eigen::Vector3d &position,
                             const long &count) {
        std::priority_queue<std::pair<double, long>> q;
        for (long i = 0; i < electronPositions.numberOfEntities(); i++) {
            q.push(std::pair<double, long>(-Metrics::distance(electronPositions[i], position),
                                           i));
        };

        std::list<long> excludedIndices = getNonValenceIndices(electronPositions, nuclei);

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
    getNearestValenceIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei, const int &index,
                             const long &count) {
        return getNearestValenceIndices(electronPositions, nuclei, nuclei[index].position(),
                                        count);
    };

    std::list<long>
    getNearestValenceIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei, const int &index1,
                             const int &index2,
                             const long &count) {
        Eigen::Vector3d position =
                (nuclei[index1].position() + nuclei[index2].position()) / 2;
        return getNearestValenceIndices(electronPositions, nuclei, position, count);
    };

    std::list<long> getNearestElectronsIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei,
                                               const Eigen::Vector3d &position,
                                               const long &count) {
        std::priority_queue<std::pair<double, int>> q;
        for (long i = 0; i < electronPositions.numberOfEntities(); i++) {
            q.push(std::pair<double, long>(-Metrics::distance(electronPositions[i], position),
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
    getNearestElectronsIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei, const int &index,
                               const long &count) {
        return getNearestElectronsIndices(electronPositions, nuclei, nuclei[index].position(),
                                          count);
    };

    std::list<long>
    getNearestElectronsIndices(const PositionsVector &electronPositions, const AtomsVector &nuclei, const int &index1,
                               const int &index2,
                               const long &count) {
        Eigen::Vector3d position =
                (nuclei[index1].position() + nuclei[index2].position()) / 2;
        return getNearestElectronsIndices(electronPositions, nuclei, position, count);
    };
}
