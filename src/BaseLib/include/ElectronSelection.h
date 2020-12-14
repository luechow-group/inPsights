// Copyright (C) 2019 Leonard Reuter.
// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ELECTRONSELECTION_H
#define INPSIGHTS_ELECTRONSELECTION_H

#include <ParticlesVector.h>
#include <ISettings.h>
#include <Metrics.h>

#include <functional>

namespace Settings {
    class ElectronSelection : public ISettings {
    public:
        Property<double> maximalDistance = {10.0, VARNAME(maximalDistance)};
        Property<long> maximalCount = {0, VARNAME(maximalCount)};
        Property<std::string> distanceMode = {"minimum", VARNAME(distanceMode)}; //TODO
        Property<bool> invertSelection = {false, VARNAME(invertSelection)};
        Property<bool> valenceOnly = {true, VARNAME(valenceOnly)};

        std::vector<Eigen::Vector3d> positions;

        ElectronSelection();
        explicit ElectronSelection(const YAML::Node &node, const AtomsVector& atoms);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::ElectronSelection)

namespace ElectronSelection {
    inline Settings::ElectronSelection settings {};

    std::list<long>
    getNonValenceIndices(const ElectronsVector &electrons, const Atom &nucleus);

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei);

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const std::vector<Eigen::Vector3d> &positions, long maximalCount,
                               const bool &valenceOnly, double maximalDistance,
                               std::function<double(const Eigen::Vector3d &,
                                                    const std::vector<Eigen::Vector3d> &)> &distanceFunction);

    std::list<long>
    getNearestPositionIndices(const PositionsVector &positions, const Eigen::Vector3d &position, long count);


    std::list<long> invertedIndices(const std::list<long>& indices, std::size_t size);

    // convenience wrapper function for clusterers
    std::list<long> getRelevantIndices(const ElectronsVector &electrons, const AtomsVector& nuclei);
}

namespace YAML{
    Eigen::Vector3d decodePosition(const YAML::Node &node, const AtomsVector &nuclei);
}

#endif //INPSIGHTS_ELECTRONSELECTION_H
