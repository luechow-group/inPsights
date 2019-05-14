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

    std::list<long> getNonValenceIndices(const MolecularGeometry &molecularGeometry, const int &k) {
        const Elements::ElementType &elementType = molecularGeometry.atoms()[k].type();
        const long coreElectronsNumber =
                Elements::ElementInfo::Z(elementType) - Elements::ElementInfo::valElectrons(elementType);
        return getNearestElectronsIndices(molecularGeometry, k, coreElectronsNumber);
    };

    std::list<long> getNonValenceIndices(const MolecularGeometry &molecularGeometry) {
        std::list<long> indices;

        for (int k = 0; k < molecularGeometry.atoms().numberOfEntities(); ++k)
            indices.splice(indices.end(), getNonValenceIndices(molecularGeometry, k));
        indices.sort();
        return indices;
    };

    Eigen::VectorXd
    getNearestValencePositions(const MolecularGeometry &molecularGeometry, const Eigen::Vector3d &position,
                               const long &count) {
        return longIndexListToPositionArrayXi(getNearestValenceIndices(molecularGeometry, position, count)).unaryExpr(
                molecularGeometry.electrons().positionsVector().asEigenVector());
    };

    Eigen::VectorXd
    getNearestValencePositions(const MolecularGeometry &molecularGeometry, const int &index, const long &count) {
        return getNearestValencePositions(molecularGeometry, molecularGeometry.atoms()[index].position(),
                                          count);
    };

    Eigen::VectorXd
    getNearestValencePositions(const MolecularGeometry &molecularGeometry, const int &index1, const int &index2,
                               const long &count) {
        Eigen::Vector3d position =
                (molecularGeometry.atoms()[index1].position() + molecularGeometry.atoms()[index2].position()) / 2;
        return getNearestValencePositions(molecularGeometry, position, count);
    };

    std::list<long>
    getNearestValenceIndices(const MolecularGeometry &molecularGeometry, const Eigen::Vector3d &position,
                             const long &count) {
        std::priority_queue<std::pair<double, long>> q;
        for (long i = 0; i < molecularGeometry.electrons().numberOfEntities(); i++) {
            q.push(std::pair<double, long>(-Metrics::distance(molecularGeometry.electrons()[i].position(), position),
                                           i));
        };

        std::list<long> excludedIndices = getNonValenceIndices(molecularGeometry);

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
    getNearestValenceIndices(const MolecularGeometry &molecularGeometry, const int &index,
                             const long &count) {
        return getNearestValenceIndices(molecularGeometry, molecularGeometry.atoms()[index].position(),
                                        count);
    };

    std::list<long>
    getNearestValenceIndices(const MolecularGeometry &molecularGeometry, const int &index1,
                             const int &index2,
                             const long &count) {
        Eigen::Vector3d position =
                (molecularGeometry.atoms()[index1].position() + molecularGeometry.atoms()[index2].position()) / 2;
        return getNearestValenceIndices(molecularGeometry, position, count);
    };

    std::list<long> getNearestElectronsIndices(const MolecularGeometry &molecularGeometry,
                                               const Eigen::Vector3d &position,
                                               const long &count) {
        std::priority_queue<std::pair<double, int>> q;
        for (long i = 0; i < molecularGeometry.electrons().numberOfEntities(); i++) {
            q.push(std::pair<double, long>(-Metrics::distance(molecularGeometry.electrons()[i].position(), position),
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
    getNearestElectronsIndices(const MolecularGeometry &molecularGeometry, const int &index,
                               const long &count) {
        return getNearestElectronsIndices(molecularGeometry, molecularGeometry.atoms()[index].position(),
                                          count);
    };

    std::list<long>
    getNearestElectronsIndices(const MolecularGeometry &molecularGeometry, const int &index1,
                               const int &index2,
                               const long &count) {
        Eigen::Vector3d position =
                (molecularGeometry.atoms()[index1].position() + molecularGeometry.atoms()[index2].position()) / 2;
        return getNearestElectronsIndices(molecularGeometry, position, count);
    };
}
