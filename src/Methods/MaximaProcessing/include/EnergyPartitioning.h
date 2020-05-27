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

#ifndef INPSIGHTS_ENERGYPARTITIONING_H
#define INPSIGHTS_ENERGYPARTITIONING_H

#include <ParticlesVector.h>
#include <MolecularGeometry.h>
#include <EnergyStatistics.h>
#include <GraphAnalysis.h>
#include <Motifs.h>
#include <Eigen/Core>

namespace EnergyPartitioning {

    double calculateTotalEnergy(const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);

    namespace MotifBased {
        double calculateSelfInteractionEnergy(const Motif &motif,
                                              const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                              const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);

        double caclulateInteractionEnergy(const Motif &motif, const Motif &otherMotif,
                                          const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                          const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);

        std::pair<Eigen::VectorXd, Eigen::MatrixXd> calculateInteractionEnergies(
                const Motifs& motifs,
                                                  const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                                  const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);
    }

    namespace ParticleBased {
        Eigen::VectorXd oneAtomEnergies(const Eigen::MatrixXd &Vnn);

        Eigen::VectorXd oneElectronEnergies(
                const Eigen::VectorXd &Te,
                const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven,
                const Eigen::MatrixXd &Vnn);
    }
}

#endif //INPSIGHTS_ENERGYPARTITIONING_H
