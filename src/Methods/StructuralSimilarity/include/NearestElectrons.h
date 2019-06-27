//
// Created by leonard on 13.05.19.
//

#ifndef INPSIGHTS_NEARESTELECTRONS_H
#define INPSIGHTS_NEARESTELECTRONS_H

#include <ParticlesVector.h>

namespace NearestElectrons {
    std::list<long>
    getNonValenceIndices(const ElectronsVector &electrons, const Atom &nucleus);

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei);

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                             const std::vector<Eigen::Vector3d> &positions,
                             const long &count);

    std::list<long>
    getNearestValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                             const Eigen::Vector3d &position,
                             const long &count);

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const Eigen::Vector3d &position,
                               const long &count);
}

#endif //INPSIGHTS_NEARESTELECTRONS_H
