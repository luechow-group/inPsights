// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "TotalWeightDifferenceAnalyzer.h"
#include "GraphAnalysis.h"
#include <ClusterAnalysis.h>
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

void TotalWeightDifferenceAnalyzer::analyze(const Cluster& cluster) {
    unsigned nIncrements = settings.increments();
    double h = settings.radiusIncrement();
    auto startRadius = settings.startRadius();

    // prepare data needed for the analysis
    std::size_t totalNumberOfMaxima = cluster.numberOfLeaves();
    auto mat = ClusterAnalysis::calculateBestMatchDistanceMatrix(cluster);
    totalWeightDifferences_.clear();


    // get weight of the clusters in cluster and put it in prev weights
    std::vector<std::set<Eigen::Index>> clusterIdsListsOfPreviousRadiusClustering;
    std::vector<double> weightsOfPreviousRadiusClustering;
    for (const auto & [clusterClusterId, cluster] : enumerate(cluster)) {
        clusterIdsListsOfPreviousRadiusClustering.emplace_back(std::set<Eigen::Index>({static_cast<Eigen::Index>(clusterClusterId)}));
        weightsOfPreviousRadiusClustering.emplace_back(static_cast<double>(cluster.numberOfLeaves()) / static_cast<double>(totalNumberOfMaxima));
    }

    for (unsigned i = 0; i <= nIncrements; ++i) {
        double radius = startRadius + i*h;

        auto adjacencyMatrix = GraphAnalysis::lowerOrEqualFilter(mat, radius);
        auto clusterIdsListsOfCurrentRadiusClustering = GraphAnalysis::findGraphClusters(adjacencyMatrix);
        auto prevToCurrClusterIdMap = GraphAnalysis::findMergeMap(clusterIdsListsOfPreviousRadiusClustering, clusterIdsListsOfCurrentRadiusClustering);
        assert(prevToCurrClusterIdMap.size() == clusterIdsListsOfPreviousRadiusClustering.size());

        double totalWeightDifference = 0.0;
        std::vector<double> weightsOfCurrentRadiusClustering;

        // iterate over clusters found for the specified radius
        for (const auto [currentClusterId, clusterIdsListOfCurrentCluster] : enumerate(clusterIdsListsOfCurrentRadiusClustering)) {

            // calculate weight of a cluster of connected cluster ids
            double currentClusterWeight = 0;
            for(auto clusterId : clusterIdsListOfCurrentCluster)
                currentClusterWeight += double(cluster[clusterId].numberOfLeaves());

            currentClusterWeight /= double(totalNumberOfMaxima);
            weightsOfCurrentRadiusClustering.emplace_back(currentClusterWeight);

            // find all clusterIds merged into the current clusterId in the map
            auto clusterIdsOfPreviousClustersMergedIntoCurrentCluster
            = MapUtils::findByValue(prevToCurrClusterIdMap, std::size_t(currentClusterId));

            assert(clusterIdsOfPreviousClustersMergedIntoCurrentCluster.size() > 0
            && "Every current cluster should have a merge from a previous one.");
            assert(clusterIdsOfPreviousClustersMergedIntoCurrentCluster.size() <= clusterIdsListOfCurrentCluster.size()
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
        clusterIdsListsOfPreviousRadiusClustering = clusterIdsListsOfCurrentRadiusClustering;
    }
}

std::vector<double> TotalWeightDifferenceAnalyzer::getResults(){
    return totalWeightDifferences_;
}
