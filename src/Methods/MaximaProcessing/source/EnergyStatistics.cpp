// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <EnergyStatistics.h>

namespace EnergyStatistics {
    ElectronicEnergy::ElectronicEnergy(
            const VectorStatistics &Te,
            const TriangularMatrixStatistics &Vee,
            const MatrixStatistics &Ven)
            :
            numberOfElectrons_(Te.mean().size()),
            Te_(Te), Vee_(Vee), Ven_(Ven) {}

    const VectorStatistics &ElectronicEnergy::Te() const {
        return Te_;
    }

    const TriangularMatrixStatistics &ElectronicEnergy::Vee() const {
        return Vee_;
    }

    const MatrixStatistics &ElectronicEnergy::Ven() const {
        return Ven_;
    }

    Eigen::Index ElectronicEnergy::numberOfElectrons() const {
        return numberOfElectrons_;
    }
}
