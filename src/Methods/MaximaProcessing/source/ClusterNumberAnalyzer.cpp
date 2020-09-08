// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <ClusterNumberAnalyzer.h>
#include <GraphAnalysis.h>
#include <ClusterAnalysis.h>
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
void  ClusterNumberAnalyzer::analyze(const Cluster& cluster) {
    unsigned nIncrements = settings.increments();
    double h = settings.radiusIncrement();
    auto startRadius = settings.startRadius();
    auto minimalWeight = settings.minimalWeight();

    auto mat = ClusterAnalysis::calculateBestMatchDistanceMatrix(cluster);
    clusterNumbers_.clear();

    std::size_t totalNumberOfMaxima = cluster.numberOfLeaves();

    std::size_t significantClusterCount = 0;

    for(const auto& subcluster : cluster){
        auto weight = static_cast<double>(subcluster.numberOfLeaves()) / static_cast<double>(totalNumberOfMaxima);
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
        for(const auto& graphCluster : graphClusters){

            // calculate weight for cluster
            std::size_t summedMaxima = 0;
            for(auto index : graphCluster)
                summedMaxima += cluster[index].numberOfLeaves();
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

