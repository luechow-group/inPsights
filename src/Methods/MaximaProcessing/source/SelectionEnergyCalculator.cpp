/* Copyright (C) 2020 Michael Heuer.
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

#include "SelectionEnergyCalculator.h"
#include <Reference.h>
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
