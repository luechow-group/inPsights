// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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

bool MolecularSelection::operator==(const MolecularSelection &rhs){
    return this->electrons_ == rhs.electrons_ && this->nuclei_ == rhs.nuclei_;
}

namespace YAML {
    Emitter &operator<<(Emitter &out, const MolecularSelection &rhs) {
        out << BeginMap
            << Key << "Nuclei" << Value << rhs.nuclei_
            << Key << "Electrons" << Value << rhs.electrons_
            << EndMap;
        return out;
    }
}

DynamicMolecularSelection::DynamicMolecularSelection(
        Settings::ParticleSelection particleSelectionSettings, ParticleIndices nucleiIndices)
        :
        particleSelectionSettings_(particleSelectionSettings),
        nuclei_(nucleiIndices) {}

ParticleIndices
DynamicMolecularSelection::selectNearestElectrons(const ElectronsVector &electrons, const AtomsVector &nuclei) const {

    std::function<double(const Eigen::Vector3d &,const std::vector<Eigen::Vector3d> &)> minimalDistanceFunction =
            Metrics::minimalDistance<2>;

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
DynamicMolecularSelection::toMolecularSelection(const ElectronsVector &electrons, const AtomsVector &nuclei) const {

    auto electronInidces = selectNearestElectrons(electrons, nuclei);

    return {electronInidces, nuclei_};
}