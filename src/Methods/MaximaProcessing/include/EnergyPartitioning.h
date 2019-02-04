//
// Created by Michael Heuer on 2019-02-03.
//

#ifndef INPSIGHTS_ENERGYPARTITIONING_H
#define INPSIGHTS_ENERGYPARTITIONING_H

#include <ParticlesVector.h>
#include <MolecularGeometry.h>
#include <EnergyStatistics.h>
#include <GraphAnalysis.h>
#include <MotifAnalysis.h>
#include <Eigen/Core>

class MotifEnergies{
public:
    void addPair(const std::pair<MotifAnalysis::Motif, MotifAnalysis::Motif>& pair, double energy){
        interactionEnergies_.emplace(std::make_pair(pair, energy));
    }

    const std::map<std::pair<MotifAnalysis::Motif, MotifAnalysis::Motif>, double>& interactionEnergies(){
        return interactionEnergies_;
    }

private:
    std::map<std::pair<MotifAnalysis::Motif, MotifAnalysis::Motif>, double> interactionEnergies_;
};


namespace EnergyPartitioning {

    // Better Use
    namespace MotifBased {

        double calculateSelfInteractionEnergy(const MotifAnalysis::Motif &motif,
                                              const EnergyStatistics::ElectronicEnergy &electronicEnergy);

        double caclulateInteractionEnergy(const MotifAnalysis::Motif &motif, const MotifAnalysis::Motif &otherMotif,
                                          const EnergyStatistics::ElectronicEnergy &electronicEnergy);

        MotifEnergies calculateInterationEnergies(const MotifAnalysis::Motifs& motif,const EnergyStatistics::ElectronicEnergy& electronicEnergy);

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
