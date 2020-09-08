// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
