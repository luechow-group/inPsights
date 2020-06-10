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

#include "TotalWeightDifferenceAnalyzer.h"
#include <GraphAnalysis.h>
#include <GroupAnalysis.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <MapUtils.h>

namespace Settings {
    TotalWeightDifferenceAnalyzer::TotalWeightDifferenceAnalyzer()
            : ISettings(VARNAME(TotalWeightDifferenceAnalyzer)) {}

    TotalWeightDifferenceAnalyzer::TotalWeightDifferenceAnalyzer(const YAML::Node &node)
            : TotalWeightDifferenceAnalyzer() {
        doubleProperty::decode(node, startRadius);
        unsignedProperty::decode(node, increments);
        doubleProperty::decode(node, radiusIncrement);
    }

    void TotalWeightDifferenceAnalyzer::appendToNode(YAML::Node &node) const {
        node[className][startRadius.name()] = startRadius();
        node[className][increments.name()] = increments();
        node[className][radiusIncrement.name()] = radiusIncrement();
    }
}
YAML_SETTINGS_DEFINITION(Settings::TotalWeightDifferenceAnalyzer)

Settings::TotalWeightDifferenceAnalyzer TotalWeightDifferenceAnalyzer::settings = Settings::TotalWeightDifferenceAnalyzer();

void TotalWeightDifferenceAnalyzer::analyze(const Group& group) {
    unsigned nIncrements = settings.increments();
    double h = settings.radiusIncrement();
    auto startRadius = settings.startRadius();

    // prepare data needed for the analysis
    std::size_t totalNumberOfMaxima = group.numberOfLeaves();
    auto mat = GroupAnalysis::calculateBestMatchDistanceMatrix(group);
    totalWeightDifferences_.clear();


    // get weight of the clusters in group and put it in prev weights
    std::vector<std::set<Eigen::Index>> groupIdsListsOfPreviousRadiusClustering;
    std::vector<double> weightsOfPreviousRadiusClustering;
    for (const auto & [clusterGroupId, cluster] : enumerate(group)) {
        groupIdsListsOfPreviousRadiusClustering.emplace_back(std::set<Eigen::Index>({static_cast<Eigen::Index>(clusterGroupId)}));
        weightsOfPreviousRadiusClustering.emplace_back(static_cast<double>(cluster.numberOfLeaves()) / static_cast<double>(totalNumberOfMaxima));
    }

    for (unsigned i = 0; i <= nIncrements; ++i) {
        double radius = startRadius + i*h;

        auto adjacencyMatrix = GraphAnalysis::lowerOrEqualFilter(mat, radius);
        auto groupIdsListsOfCurrentRadiusClustering = GraphAnalysis::findGraphClusters(adjacencyMatrix);
        auto prevToCurrClusterIdMap = GraphAnalysis::findMergeMap(groupIdsListsOfPreviousRadiusClustering, groupIdsListsOfCurrentRadiusClustering);
        assert(prevToCurrClusterIdMap.size() == groupIdsListsOfPreviousRadiusClustering.size());

        double totalWeightDifference = 0.0;
        std::vector<double> weightsOfCurrentRadiusClustering;

        // iterate over clusters found for the specified radius
        for (const auto [currentClusterId, groupIdsListOfCurrentCluster] : enumerate(groupIdsListsOfCurrentRadiusClustering)) {

            // calculate weight of a cluster of connected group ids
            double currentClusterWeight = 0;
            for(auto groupId : groupIdsListOfCurrentCluster)
                currentClusterWeight += double(group[groupId].numberOfLeaves());

            currentClusterWeight /= double(totalNumberOfMaxima);
            weightsOfCurrentRadiusClustering.emplace_back(currentClusterWeight);

            // find all clusterIds merged into the current clusterId in the map
            auto clusterIdsOfPreviousClustersMergedIntoCurrentCluster
            = MapUtils::findByValue(prevToCurrClusterIdMap, std::size_t(currentClusterId));

            assert(clusterIdsOfPreviousClustersMergedIntoCurrentCluster.size() > 0
            && "Every current cluster should have a merge from a previous one.");
            assert(clusterIdsOfPreviousClustersMergedIntoCurrentCluster.size() <= groupIdsListOfCurrentCluster.size()
            && "Maximally, all previous clusters can be merged into the current one.");

            // find clusterId of previous cluster with the largest weight
            auto heaviestPreviousClusterId = std::max_element(
                    clusterIdsOfPreviousClustersMergedIntoCurrentCluster.begin(),
                    clusterIdsOfPreviousClustersMergedIntoCurrentCluster.end(),
                    [weightsOfPreviousRadiusClustering] (std::size_t lhs, size_t rhs) {
                        return weightsOfPreviousRadiusClustering[lhs] < weightsOfPreviousRadiusClustering[rhs];
            });
            
            // remove it
            clusterIdsOfPreviousClustersMergedIntoCurrentCluster.erase(heaviestPreviousClusterId);

            // calculate weight difference
            double weightMergedIntoCurrentCluster = 0.0;
            for(auto previousClusterIdMergedIntoCurrent : clusterIdsOfPreviousClustersMergedIntoCurrentCluster) {
                weightMergedIntoCurrentCluster += weightsOfPreviousRadiusClustering[previousClusterIdMergedIntoCurrent];
            }

            totalWeightDifference += weightMergedIntoCurrentCluster;

        }
        totalWeightDifferences_.emplace_back(totalWeightDifference);
        weightsOfPreviousRadiusClustering = weightsOfCurrentRadiusClustering;
        groupIdsListsOfPreviousRadiusClustering = groupIdsListsOfCurrentRadiusClustering;
    }
}

std::vector<double> TotalWeightDifferenceAnalyzer::getResults(){
    return totalWeightDifferences_;
}
