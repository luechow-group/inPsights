//
// Created by Michael Heuer on 2019-02-03.
//

#ifndef INPSIGHTS_ENERGYPARTITIONING_H
#define INPSIGHTS_ENERGYPARTITIONING_H

#include <ParticlesVector.h>
#include <MolecularGeometry.h>
#include <EnergyStatistics.h>
#include <GraphAnalysis.h>
#include <MotifEnergies.h>
#include <Motifs.h>
#include <Eigen/Core>

namespace EnergyPartitioning {

    // Better Use
    namespace MotifBased {

        double calculateSelfInteractionEnergy(const Motif &motif,
                                              const EnergyStatistics::ElectronicEnergy &electronicEnergy);

        double caclulateInteractionEnergy(const Motif &motif, const Motif &otherMotif,
                                          const EnergyStatistics::ElectronicEnergy &electronicEnergy);

        MotifEnergies calculateInterationEnergies(const Motifs& motif,const EnergyStatistics::ElectronicEnergy& electronicEnergy);

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
