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

#include <ClusterNumberAnalyzer.h>
#include <GraphAnalysis.h>
#include <GroupAnalysis.h>
#include <ValueSorter.h>
#include <spdlog/spdlog.h>
#include <ErrorHandling.h>

namespace Settings {
    ClusterNumberAnalyzer::ClusterNumberAnalyzer()
            : ISettings(VARNAME(ClusterNumberAnalyzer)) {}

    ClusterNumberAnalyzer::ClusterNumberAnalyzer(const YAML::Node &node)
            : ClusterNumberAnalyzer() {
        doubleProperty::decode(node, startRadius);
        doubleProperty::decode(node, radiusIncrement);
        unsignedProperty::decode(node, increments);
        doubleProperty::decode(node, minimalWeight);
    }

    void ClusterNumberAnalyzer::appendToNode(YAML::Node &node) const {
        node[className][startRadius.name()] = startRadius();
        node[className][radiusIncrement.name()] = radiusIncrement();
        node[className][increments.name()] = increments();
        node[className][minimalWeight.name()] = minimalWeight();
    }
}
YAML_SETTINGS_DEFINITION(Settings::ClusterNumberAnalyzer)

Settings::ClusterNumberAnalyzer ClusterNumberAnalyzer::settings = Settings::ClusterNumberAnalyzer();

#include <iomanip>
void  ClusterNumberAnalyzer::analyze(const Group& group) {
    unsigned nIncrements = settings.increments();
    double h = settings.radiusIncrement();
    auto startRadius = settings.startRadius();
    auto minimalWeight = settings.minimalWeight();

    auto mat = GroupAnalysis::calculateBestMatchDistanceMatrix(group);
    clusterNumbers_.clear();

    std::size_t totalNumberOfMaxima = group.numberOfLeaves();

    std::size_t significantClusterCount = 0;

    for(const auto& subgroup : group){
        auto weight = static_cast<double>(subgroup.numberOfLeaves()) / static_cast<double>(totalNumberOfMaxima);
        if(weight > minimalWeight)
            significantClusterCount++;
    }
    clusterNumbers_.emplace_back(significantClusterCount);

    for (unsigned i = 0; i <= nIncrements; ++i) {
        double radius = startRadius + i*h;
        auto adjacencyMatrix = GraphAnalysis::lowerOrEqualFilter(mat, radius);
        auto graphClusters = GraphAnalysis::findGraphClusters(adjacencyMatrix);

        // filter the number of significant clusters in the set of found clusters
        significantClusterCount = 0;
        for(const auto& cluster : graphClusters){

            // calculate weight for cluster
            std::size_t summedMaxima = 0;
            for(auto index : cluster)
                summedMaxima += group[index].numberOfLeaves();
            auto weight = static_cast<double>(summedMaxima) / static_cast<double>(totalNumberOfMaxima);

            // check if weight is significant enough
            if(weight > minimalWeight)
                significantClusterCount++;
        }
        clusterNumbers_.emplace_back(significantClusterCount);
    }
}

std::vector<std::size_t> ClusterNumberAnalyzer::getResults() {
    return clusterNumbers_;
}

