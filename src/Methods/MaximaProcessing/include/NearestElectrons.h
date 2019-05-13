//
// Created by leonard on 13.05.19.
//

#ifndef INPSIGHTS_NEARESTELECTRONS_H
#define INPSIGHTS_NEARESTELECTRONS_H

#include <MolecularGeometry.h>

namespace NearestElectrons {

    std::list<long> getNonValenceIndices(const MolecularGeometry &molecularGeometry, const int &k);

    std::list<long> getNonValenceIndices(const MolecularGeometry &molecularGeometry);

    std::list<long>
    getNearestValenceIndices(const MolecularGeometry &molecularGeometry, const Eigen::Vector3d &position,
                             const long &count);

    std::list<long>
    getNearestValenceIndices(const MolecularGeometry &molecularGeometry, const int &index, const long &count);

    std::list<long>
    getNearestValenceIndices(const MolecularGeometry &molecularGeometry, const int &index1, const int &index2,
                             const long &count);

    std::list<long>
    getNearestElectronsIndices(const MolecularGeometry &molecularGeometry, const Eigen::Vector3d &position,
                               const long &count);

    std::list<long>
    getNearestElectronsIndices(const MolecularGeometry &molecularGeometry, const int &index, const long &count);

    std::list<long>
    getNearestElectronsIndices(const MolecularGeometry &molecularGeometry, const int &index1, const int &index2,
                               const long &count);
}

#endif //INPSIGHTS_NEARESTELECTRONS_H
