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

#include <MolecularSelection.h>

#include <utility>


MolecularSelection::MolecularSelection(ParticleIndices electronIndices)
    :
    electrons_(std::move(electronIndices)),
    nuclei_()
    {}


MolecularSelection::MolecularSelection(ParticleIndices electronIndices, ParticleIndices nucleiIndices)
    :
    electrons_(std::move(electronIndices)),
    nuclei_(std::move(nucleiIndices))
    {}

DynamicMolecularSelection::DynamicMolecularSelection(
        Settings::ParticleSelection particleSelectionSettings, ParticleIndices nucleiIndices)
        :
        particleSelectionSettings_(particleSelectionSettings),
        nuclei_(nucleiIndices) {}

ParticleIndices
DynamicMolecularSelection::selectNearestElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei) {

    std::function<double(const Eigen::Vector3d &,const std::vector<Eigen::Vector3d> &)> minimalDistanceFunction = Metrics::minimalDistance<2>;

    auto selectedElectrons = ParticleSelection::getNearestElectronsIndices(electrons, nuclei,
                                                                           particleSelectionSettings_.positions,
                                                                           particleSelectionSettings_.maximalCount(),
                                                                           particleSelectionSettings_.valenceOnly(),
                                                                           particleSelectionSettings_.maximalDistance(),
                                                                           minimalDistanceFunction
    );

    ParticleIndices indices;

    for(auto i : selectedElectrons)
        indices.indices().emplace(i);

    return indices;
}

MolecularSelection
DynamicMolecularSelection::toMolecularSelection(const ElectronsVector &electrons, const AtomsVector &nuclei) {

    auto electronInidces = selectNearestElectrons(electrons, nuclei);

    return {electronInidces, nuclei_};
}