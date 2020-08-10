/* Copyright (C) 2019 Leonard Reuter.
 * Copyright (C) 2020 Michael Heuer.
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

#ifndef INPSIGHTS_PARTICLESELECTION_H
#define INPSIGHTS_PARTICLESELECTION_H

#include <ParticlesVector.h>
#include <ISettings.h>
#include <Metrics.h>

#include <functional>

namespace Settings {
    class ParticleSelection : public ISettings {
    public:
        Property<double> maximalDistance = {10.0, VARNAME(maximalDistance)};
        Property<long> maximalCount = {0, VARNAME(maximalCount)};
        Property<std::string> distanceMode = {"minimum", VARNAME(distanceMode)}; //TODO
        Property<bool> invertSelection = {false, VARNAME(invertSelection)};
        Property<bool> valenceOnly = {true, VARNAME(valenceOnly)};

        AtomsVector atoms;
        std::vector<Eigen::Vector3d> positions;

        ParticleSelection(const AtomsVector& atoms = {});
        explicit ParticleSelection(const YAML::Node &node, const AtomsVector& atoms = {});
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::ParticleSelection)

namespace ParticleSelection {
    inline Settings::ParticleSelection settings {};

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
    std::list<long> getRelevantIndices(const ElectronsVector &electrons);
}

#endif //INPSIGHTS_PARTICLESELECTION_H
