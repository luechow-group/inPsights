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

#include <GraphClusterer.h>

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

#include <GraphClusterer.h>
#include <GraphAnalysis.h>
#include <BestMatchDistance.h>
#include <ValueSorter.h>
#include <spdlog/spdlog.h>
#include <ErrorHandling.h>

namespace Settings {
    GraphClusterer::GraphClusterer()
            : ISettings(VARNAME(GraphClusterer)) {}

    GraphClusterer::GraphClusterer(const YAML::Node &node)
            : GraphClusterer() {
        doubleProperty::decode(node, startRadius);
        doubleProperty::decode(node, endRadius);
        doubleProperty::decode(node, radiusIncrement);
        doubleProperty::decode(node, minimalWeight);
    }

    void GraphClusterer::appendToNode(YAML::Node &node) const {
        node[className][startRadius.name()] = startRadius();
        node[className][endRadius.name()] = endRadius();
        node[className][radiusIncrement.name()] = radiusIncrement();
        node[className][minimalWeight.name()] = minimalWeight();
    }
}
YAML_SETTINGS_DEFINITION(Settings::GraphClusterer)

Settings::GraphClusterer GraphClusterer::settings = Settings::GraphClusterer();


GraphClusterer::GraphClusterer(Group& group)
        : mat_(calculateAdjacencyMatrix(group)) {}

// constructs an adjacency matrix with indices being determined from the order in the group
Eigen::MatrixXd  GraphClusterer::calculateAdjacencyMatrix(Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    // first, sort references by value
    group.sortAll();

    // construct matrix
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(group.size(),group.size());
    Eigen::Index nClusters = group.size();
    spdlog::info("number of clusters {}", nClusters);
    for (Eigen::Index i = 0; i < nClusters-1; ++i) {
        for (Eigen::Index j = i+1; j < nClusters; ++j) {

            auto [norm, perm] = BestMatch::Distance::compare<Spin, Eigen::Infinity, 2>(
                    group[i].representative()->maximum(),
                    group[j].representative()->maximum());
            mat(i,j) = norm;
        }
    }
    // symmetrization
    return mat.selfadjointView<Eigen::Upper>();
}

void GraphClusterer::cluster(Group& group) {
    throw NotImplemented();
}

std::vector<std::size_t >  GraphClusterer::scanClusterSizeWithDistance(const Group& group) {
    double b = settings.endRadius();
    double h = settings.radiusIncrement();
    auto dist = settings.startRadius();
    auto minimalWeight = settings.minimalWeight();

    std::size_t totalNumberOfMaxima = group.numberOfLeaves();

    std::vector<std::size_t> clusterSizes;
    std::size_t significantClusterCount = 0;

    for(const auto& subgroup : group){
        auto weight = static_cast<double>(subgroup.numberOfLeaves()) / static_cast<double>(totalNumberOfMaxima);
        if(weight > minimalWeight)
            significantClusterCount++;
    }
    clusterSizes.emplace_back(significantClusterCount);

    while( dist <= b){
        dist += h;
        auto adjacencyMatrix = GraphAnalysis::lowerOrEqualFilter(mat_, dist);
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
        clusterSizes.emplace_back(significantClusterCount);


    }
    return clusterSizes;
}
