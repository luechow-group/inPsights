/* Copyright 2020 heuer
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

/* Copyright 2020 heuer
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

/* Copyright (C) 2018-2019 Michael Heuer.
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
    auto mat = GroupAnalysis::calculateAdjacencyMatrix(group);
    totalWeightDifferences_.clear();


    // get weight of the clusters in group and put it in prev weights
    std::vector<std::list<Eigen::Index>> previousClusters;
    std::vector<double> prevWeights;
    for (const auto & [clusterIndex, cluster] : enumerate(group)) {
        previousClusters.emplace_back(std::list<Eigen::Index>({static_cast<Eigen::Index>(clusterIndex)}));
        prevWeights.emplace_back(static_cast<double>(cluster.numberOfLeaves()) / static_cast<double>(totalNumberOfMaxima));
    }

    for (unsigned i = 0; i <= nIncrements; ++i) {
        double radius = startRadius + i*h;

        auto adjacencyMatrix = GraphAnalysis::lowerOrEqualFilter(mat, radius);
        auto currentClustersMembersIndices = GraphAnalysis::findGraphClusters(adjacencyMatrix);
        auto prevToCurrMap = GraphAnalysis::findMergeMap(previousClusters, currentClustersMembersIndices);

        double totalWeightDifference = 0.0;
        std::vector<double> newWeights;

        for (const auto& [currentClusterIndex, currentClusterMemberIndices] : enumerate(currentClustersMembersIndices)) {

            double currentWeight = 0;
            for(auto memberIndex : currentClusterMemberIndices)
                currentWeight += group[memberIndex].numberOfLeaves();
            currentWeight /= double(totalNumberOfMaxima);
            newWeights.emplace_back(currentWeight);

            // find all prev clusters that were merged into the current cluster
            std::vector<Eigen::Index> vec;

            // find all members in the previous map currentClusterMemberIndices
            bool foundQ = GraphAnalysis::findByValue(vec, prevToCurrMap, Eigen::Index(currentClusterIndex)); //TODO°!!!!
            assert(foundQ);
            assert(vec.size() > 0 && vec.size() <= currentClusterMemberIndices.size());

            // determine max weights of those
            double maxPrevWeight = 0.0;
            for(auto prevIndex : vec) {
                auto weight = prevWeights[prevIndex];
                if(weight > maxPrevWeight)
                    maxPrevWeight = weight;
            }
            totalWeightDifference += currentWeight - maxPrevWeight;
        }

        totalWeightDifferences_.emplace_back(totalWeightDifference);
        prevWeights = newWeights;
    }
}

std::vector<double> TotalWeightDifferenceAnalyzer::getResults(){
    return totalWeightDifferences_;
}
