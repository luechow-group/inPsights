// Copyright (C) 2019 Leonard Reuter.
// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <ParticleSelection.h>
#include <Metrics.h>
#include <ElementInfo.h>
#include <queue>
#include <algorithm>

namespace Settings {
    ParticleSelection::ParticleSelection(const AtomsVector& atoms)
            :
            ISettings(VARNAME(ParticleSelection)),
            atoms(atoms){
        distanceMode.onChange_.connect(
                [&](const std::string& value) {
                    if (value != "minimum" || value != "average")
                        throw std::invalid_argument("The distanceMode has to be minimum or average."); //TODO ADD ENUM
                });
        maximalDistance.onChange_.connect(
                [&](double value) {
                    if (value <= 0.0)
                        throw std::invalid_argument("The maximalDistance has to be larger than zero.");
                });
        maximalCount.onChange_.connect(
                [&](long value) {
                    if (value <= 0)
                        throw std::invalid_argument("The maximalCount has to be larger than zero.");
                });
    };

    ParticleSelection::ParticleSelection(const YAML::Node &node, const AtomsVector& atoms)
            : ParticleSelection(atoms) {
        doubleProperty::decode(node, maximalDistance);
        longProperty::decode(node, maximalCount);
        stringProperty::decode(node, distanceMode);
        boolProperty::decode(node, valenceOnly);
        boolProperty::decode(node, invertSelection);

        auto positionNodes = node[VARNAME(positions)];
        for (const auto &positionNode : positionNodes){
            positions.emplace_back(YAML::decodePosition(positionNode, atoms));
        }
        spdlog::info("Using the following positions:");
        for (const auto &position : positions){
            spdlog::info("{} {} {}", position[0], position[1], position[2]);
        }
    };

    void ParticleSelection::appendToNode(YAML::Node &node) const {
        node[className][maximalDistance.name()] = maximalDistance();
        node[className][maximalCount.name()] = maximalCount();
        node[className][distanceMode.name()] = distanceMode();
        node[className][valenceOnly.name()] = valenceOnly();
        node[className][invertSelection.name()] = invertSelection();

        for(const auto & p : positions)
            node[className][VARNAME(positions)].push_back(p);
    };
}

namespace ParticleSelection {
    std::list<long>
    getNonValenceIndices(const ElectronsVector &electrons, const Atom &nucleus) {
        // returns the core electron indices of a given nucleus
        const Elements::ElementType &elementType = nucleus.type();
        const long coreElectronsNumber =
                Elements::ElementInfo::Z(elementType) - Elements::ElementInfo::valenceElectrons(elementType);
        return getNearestPositionIndices(electrons.positionsVector(), nucleus.position(), coreElectronsNumber);
    };

    std::list<long> getNonValenceIndices(const ElectronsVector &electrons, const AtomsVector &nuclei) {
        // returns all core electron indices of a given vector of nuclei
        std::list<long> indices;

        for (long k = 0; k < nuclei.numberOfEntities(); ++k)
            indices.splice(indices.end(), getNonValenceIndices(electrons, nuclei[k]));
        indices.sort();
        return indices;
    };

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const std::vector<Eigen::Vector3d> &positions, long maximalCount,
                               const bool &valenceOnly, double maximalDistance,
                               std::function<double(const Eigen::Vector3d &,
                                                    const std::vector<Eigen::Vector3d> &)> &distanceFunction) {
        // returns indices of electrons closest to a vector of 'positions'
        // the metric is defined by the 'distanceFunction'
        // the indices are chosen beginning with the closest one, until 'maximalCount' or 'maximalDistance' is reached
        std::priority_queue<
                std::pair<double, int>,
                std::vector<std::pair<double, int>>,
                std::greater<std::pair<double, int>>> q;

        for (long i = 0; i < electrons.numberOfEntities(); i++) {
            q.push({distanceFunction(electrons[i].position(), positions), i});
        };

        std::list<long> excludedIndices;
        if (valenceOnly) excludedIndices = getNonValenceIndices(electrons, nuclei);

        //TODO refactor : very similar to getNearestElectronsIndices(const ElectronsVector &electrons,
        //                                                           const Eigen::Vector3d &position,
        //                                                           const long &count)
        std::list<long> indices;
        for (long j = 0; j < maximalCount; ++j) {
            // pop queue if index in pair q.top() is in excludedIndices
            while (std::find(excludedIndices.begin(), excludedIndices.end(), q.top().second) != excludedIndices.end()) {
                q.pop();
            };

            if (q.top().first > maximalDistance) break;
            indices.emplace_back(q.top().second);
            q.pop();
        };

        return indices;
    };

    std::list<long> getNearestPositionIndices(const PositionsVector &positions,
                                              const Eigen::Vector3d &position,
                                              long count) {
        // returns indices of the 'count' positions closest to 'position' (without restrictions)
        std::priority_queue<
                std::pair<double, int>,
                std::vector<std::pair<double, int>>,
                std::greater<std::pair<double, int>>> q;

        for (long i = 0; i < positions.numberOfEntities(); i++) {
            q.push({Metrics::distance(positions[i], position),i});
        };

        std::list<long> indices;
        for (int i = 0; i < count; ++i) {
            indices.emplace_back(q.top().second);
            q.pop();
        };

        indices.sort();

        return indices;
    };

    std::list<long> getRelevantIndices(const ElectronsVector &electrons) {

        auto nuclei = settings.atoms;
        auto positions = settings.positions;

        std::function<double(const Eigen::Vector3d &, const std::vector<Eigen::Vector3d> &)> distanceFunction;

        if (settings.distanceMode() == "average") { //TODO use enum
            distanceFunction = Metrics::averageDistance<2>;
        } else if (settings.distanceMode() == "minimum") {
            distanceFunction = Metrics::minimalDistance<2>;
        }

        auto indices = getNearestElectronsIndices(electrons, nuclei, positions,
                                                  settings.maximalCount(), settings.valenceOnly(),
                                                  settings.maximalDistance(),
                                                  distanceFunction);

        if (settings.invertSelection()) {
            if (settings.valenceOnly()) {
                // since selection should be inverted, core indices have to be added to 'subIndices' before inverting
                indices.splice(indices.end(), getNonValenceIndices(electrons, nuclei));
            }

            auto res = invertedIndices(indices, electrons.numberOfEntities());

            return res;
            // sorting will take the front indices. Since the invertSelection is true,
        } else
            return indices;
    }

    std::list<long> invertedIndices(const std::list<long>& indices, std::size_t size){
        assert(!indices.empty() && size > 0);
        assert(*std::min_element(indices.begin(),indices.end()) >= 0);
        assert(*std::max_element(indices.begin(),indices.end()) < long(size));

        // make sure the indices are sorted (set_difference assumes a set which is ordered by definition)
        auto sortedIndices = indices;
        sortedIndices.sort();

        std::list<long> inverted{}, all(size);
        std::iota(all.begin(), all.end(), 0);

        std::set_difference(all.begin(), all.end(), sortedIndices.begin(), sortedIndices.end(),
                            std::inserter(inverted, inverted.begin()));

        return inverted;
    };
}

namespace YAML{
    Eigen::Vector3d decodePosition(const YAML::Node &node, const AtomsVector &nuclei){
        if (node.size() != 1) throw std::invalid_argument("Nodes in positions must have size 1.");

        if (node["atAtom"].IsDefined()){
            return nuclei[node["atAtom"].as<int>()].position();
        }
        else if (node["atCoordinates"].IsDefined()) {
            return node["atCoordinates"].as<Eigen::Vector3d>();
        }
        else if (node["add"].IsDefined()) {
            Eigen::Vector3d position({0.0,0.0,0.0});
            for (const auto &positionNode : node["add"])
                position += decodePosition(positionNode, nuclei);
            return position;
        }
        else if (node["multiply"].IsDefined()) {
            Eigen::Vector3d position({0.0,0.0,0.0});
            for (const auto &positionNode : node["multiply"])
                position = position.cwiseProduct(decodePosition(positionNode, nuclei));
            return position;
        }
        else if (node["inBetween"].IsDefined()) {
            Eigen::Vector3d position({0.0,0.0,0.0});
            int positionsNumber = 0;
            for (const auto &positionNode : node["inBetween"]){
                position += decodePosition(positionNode, nuclei);
                positionsNumber++;
            }
            return position / positionsNumber;
        }
        else{
            throw std::invalid_argument("Invalid key in positions.");
        }
    }
}