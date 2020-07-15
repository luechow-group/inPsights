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

SelectionEnergyCalculator::SelectionEnergyCalculator(
        const std::vector<Sample> &samples,
        const std::vector<DynamicMolecularSelection> &selections)
        :
        samples_(samples),
        dynamicSelections_(selections),
        selectionInteractions_()
        {}

void SelectionEnergyCalculator::add(const Cluster &cluster) {

    if(cluster.isLeaf()) {
        const auto &ref = *cluster.representative();
        auto permutedNuclei = ref.nuclei();
        permutedNuclei.permute(ref.nuclearPermutation());
        auto Vnn = CoulombPotential::energies(permutedNuclei);

        std::vector<MolecularSelection> molecularSelections;
        for(auto dynamicSelection : dynamicSelections_)
            molecularSelections.emplace_back(dynamicSelection.toMolecularSelection(ref.maximum(), permutedNuclei));

        for (const auto &id : ref.sampleIds()) {
            auto &electrons = samples_[id].sample_;
            Eigen::VectorXd Te = samples_[id].kineticEnergies_;
            Eigen::MatrixXd Vee = CoulombPotential::energies(electrons);
            Eigen::MatrixXd Ven = CoulombPotential::energies(electrons, permutedNuclei);

            auto res = EnergyPartitioning::MolecularSelectionBased::calculateInteractionEnergies(
                    molecularSelections, Te, Vee, Ven, Vnn);

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
