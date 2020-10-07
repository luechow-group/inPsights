// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SelectionEnergyCalculator.h"
#include <Maximum.h>
#include <algorithm>

SelectionEnergyCalculator::SelectionEnergyCalculator(
        const std::vector<Sample> &samples,
        const std::vector<DynamicMolecularSelection> &selections)
        :
        samples_(samples),
        dynamicSelections_(selections),
        selectionInteractions_()
        {}

void SelectionEnergyCalculator::addTopLevel(const Cluster &cluster) {
    molecularSelections_ = getMolecularSelection(cluster);

    add(cluster);
}

std::vector<MolecularSelection> SelectionEnergyCalculator::getMolecularSelection(const Cluster &cluster) const {
    const auto &ref = *cluster.representative();
    auto permutedNuclei = ref.nuclei();
    permutedNuclei.permute(ref.nuclearPermutation());

    std::vector<MolecularSelection> currentMolecularSelections;
    for(const auto& dynamicSelection : this->dynamicSelections_) {
        auto sel = dynamicSelection.toMolecularSelection(ref.maximum(), permutedNuclei);

        // if specified, invert selection
        if(dynamicSelection.particleSelectionSettings_.invertSelection()){
            std::set<long> all, diff;
            generate_n(inserter(all, all.end()), ref.maximum().numberOfEntities(), [&]{ return all.size(); });

            std::set_difference(
                    all.begin(), all.end(),
                    sel.electrons_.indices().begin(), sel.electrons_.indices().end(),
                    std::inserter(diff, diff.end()));

            sel.electrons_.setIndices(diff);
        }

        currentMolecularSelections.emplace_back(sel);
    }

    return currentMolecularSelections;
}

void SelectionEnergyCalculator::add(const Cluster &cluster) {

    auto currentMolecularSelections = getMolecularSelection(cluster);
    for(const auto& [i, sel] : enumerate(molecularSelections_)) {
        assert(currentMolecularSelections[i] == sel && "Selections within a cluster should be identical.");
    }

    if(cluster.isLeaf()) {
        const auto &ref = *cluster.representative();
        auto permutedNuclei = ref.nuclei();
        permutedNuclei.permute(ref.nuclearPermutation());
        auto Vnn = CoulombPotential::energies(permutedNuclei);


        for (const auto &id : ref.sampleIds()) {
            auto &electrons = samples_[id].sample_;
            Eigen::VectorXd Te = samples_[id].kineticEnergies_;
            Eigen::MatrixXd Vee = CoulombPotential::energies(electrons);
            Eigen::MatrixXd Ven = CoulombPotential::energies(electrons, permutedNuclei);

            auto res = EnergyPartitioning::MolecularSelectionBased::calculateInteractionEnergies(
                    molecularSelections_, Te, Vee, Ven, Vnn);

            selectionInteractions_.intraEnergies.E.add(res.first.E);
            selectionInteractions_.intraEnergies.Te.add(res.first.Te);
            selectionInteractions_.intraEnergies.Vee.add(res.first.Vee);
            selectionInteractions_.intraEnergies.Ven.add(res.first.Ven);
            selectionInteractions_.intraEnergies.Vnn.add(res.first.Vnn);

            selectionInteractions_.interEnergies.E.add(res.second.E);
            selectionInteractions_.interEnergies.Vee.add(res.second.Vee);
            selectionInteractions_.interEnergies.Ven.add(res.second.Ven);
            selectionInteractions_.interEnergies.Vnn.add(res.second.Vnn);
        }
    } else {
        for (const auto& subcluster : cluster)
            add(subcluster);
    }
}

namespace YAML {
    Emitter &operator<<(Emitter &out, const SelectionEnergyCalculator::SelectionInteractionEnergies &rhs) {
        out << BeginMap
            << Key << "Intra" << Value << rhs.intraEnergies << Newline
            << Key << "Inter" << Value << rhs.interEnergies
            << EndMap;
        return out;
    }
}
