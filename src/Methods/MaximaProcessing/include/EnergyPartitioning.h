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
#include <EnergyResultsBundle.h>

namespace EnergyPartitioning {

    double calculateTotalEnergy(const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);

    namespace MolecularSelectionBased {

        std::pair<EnergyResultsBundle<Eigen::VectorXd>, EnergyResultsBundle<Eigen::MatrixXd>> calculateInteractionEnergies(
                const Motifs& motifs,
                const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);


        std::pair<EnergyResultsBundle<Eigen::VectorXd>, EnergyResultsBundle<Eigen::MatrixXd>> calculateInteractionEnergies(
                const std::vector<MolecularSelection>& selections,
                const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);

        double calculateIntraTe(const MolecularSelection &sel, const Eigen::VectorXd &Vee);

        double calculateIntraVee(const MolecularSelection &sel, const Eigen::MatrixXd &Vee);
        double calculateInterVee(const MolecularSelection &selA, const MolecularSelection &selB, const Eigen::MatrixXd &Vee);

        double calculateIntraVen(const MolecularSelection &sel, const Eigen::MatrixXd &Ven);
        double calculateInterVen(const MolecularSelection &selA, const MolecularSelection &selB, const Eigen::MatrixXd &Ven);

        double calculateIntraVnn(const MolecularSelection &sel, const Eigen::MatrixXd &Vnn);
        double calculateInterVnn(const MolecularSelection &selA, const MolecularSelection &selB, const Eigen::MatrixXd &Vnn);

        EnergyResultsBundle<double> calculateSelfInteractionEnergyBundle(
                const MolecularSelection &sel,
                const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);

        EnergyResultsBundle<double> calculateInteractionEnergyBundle(
                const MolecularSelection &selA, const MolecularSelection &selB,
                const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn);
    }

    namespace ParticleBased { // TODO remove
        Eigen::VectorXd oneAtomEnergies(const Eigen::MatrixXd &Vnn);

        Eigen::VectorXd oneElectronEnergies(
                const Eigen::VectorXd &Te,
                const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven,
                const Eigen::MatrixXd &Vnn);
    }
}

#endif //INPSIGHTS_ENERGYPARTITIONING_H
