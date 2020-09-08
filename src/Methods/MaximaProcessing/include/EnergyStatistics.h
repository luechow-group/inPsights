// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ENERGYSTATISTICS_H
#define INPSIGHTS_ENERGYSTATISTICS_H

#include "Statistics.h"

namespace EnergyStatistics {

    class ElectronicEnergy {
    public:
        ElectronicEnergy() = default;
        ElectronicEnergy(const VectorStatistics &Te,
                         const TriangularMatrixStatistics &Vee,
                         const MatrixStatistics &Ven);

        Eigen::Index numberOfElectrons() const;

        const VectorStatistics &Te() const;

        const TriangularMatrixStatistics &Vee() const;

        const MatrixStatistics &Ven() const;

    private:
        Eigen::Index numberOfElectrons_;
        VectorStatistics Te_; // kinetic energy
        TriangularMatrixStatistics Vee_; // electron-electron interaction energy
        MatrixStatistics Ven_; // electron-core interaction energy
    };


}

#endif //INPSIGHTS_ENERGYSTATISTICS_H
