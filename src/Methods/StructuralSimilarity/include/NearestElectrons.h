//
// Created by leonard on 13.05.19.
//

#ifndef INPSIGHTS_NEARESTELECTRONS_H
#define INPSIGHTS_NEARESTELECTRONS_H

#include <ParticlesVector.h>

namespace NearestElectrons {
    std::list<long>
    getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &k);

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei);

    ElectronsVector
    getNearestValenceElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const Eigen::Vector3d &position,
                               const long &count);

    ElectronsVector
    getNearestValenceElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index,
                               const long &count);

    ElectronsVector
    getNearestValenceElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index1,
                               const int &index2,
                               const long &count);

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                             const Eigen::Vector3d &position,
                             const long &count);

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index,
                             const long &count);

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index1,
                             const int &index2,
                             const long &count);

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const Eigen::Vector3d &position,
                               const long &count);

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index,
                               const long &count);

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei, const int &index1,
                               const int &index2,
                               const long &count);
}

#endif //INPSIGHTS_NEARESTELECTRONS_H
