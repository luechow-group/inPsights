/* Copyright (C) 2019 Michael Heuer.
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
