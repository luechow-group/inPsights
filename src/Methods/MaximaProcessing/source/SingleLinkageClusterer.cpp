// Copyright (C) 2018-2020 Michael Heuer.
// Copyright (C) 2018 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <SingleLinkageClusterer.h>
#include <PreClusterer.h>
#include <DistanceBasedMetric.h>
#include <Enumerate.h>
#include <Maximum.h>
#include "ClusteringMetric.h"

namespace Settings {
    SingleLinkageClusterer::SingleLinkageClusterer()
    : ISettings(VARNAME(SingleLinkageClusterer)) {
        radius.onChange_.connect(
                [&](double value) {
                    if(value < ::PreClusterer::settings.radius())
                        throw std::invalid_argument(
                                "The " + radius.name() + "=" + std::to_string(radius())
                                + " is smaller than the " + ::PreClusterer::settings.name()
                                + "::" + ::PreClusterer::settings.radius.name()
                                + " of " + std::to_string(::PreClusterer::settings.radius()) + ".");
                });
        minimalClusterSize.onChange_.connect(
                [&](size_t value) {
                    if(value < 1)
                        throw std::invalid_argument(
                                "The " + minimalClusterSize.name() + "=" + std::to_string(minimalClusterSize())
                                + " must be 1 or greater.");
                });
        local.onChange_.connect(
                [&](bool value) {
                    if (value) spdlog::info("Clustering is local");
                    else spdlog::info("Clustering is global.");
                });
    }

    SingleLinkageClusterer::SingleLinkageClusterer(const YAML::Node &node)
            : SingleLinkageClusterer() {
        doubleProperty::decode(node, radius);
        size_tProperty ::decode(node, minimalClusterSize);
        boolProperty::decode(node, local);
        boolProperty::decode(node, sortRemainder);
    }

    void SingleLinkageClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][minimalClusterSize.name()] = minimalClusterSize();
        node[className][local.name()] = local();
        node[className][sortRemainder.name()] = sortRemainder();
    }
}
YAML_SETTINGS_DEFINITION(Settings::SingleLinkageClusterer)

Settings::SingleLinkageClusterer SingleLinkageClusterer::settings = Settings::SingleLinkageClusterer();


SingleLinkageClusterer::SingleLinkageClusterer(std::vector<Sample> &samples)
        : IClusterer(samples) {};

void SingleLinkageClusterer::cluster(Cluster& cluster) {
    assert(!cluster.empty() && "The cluster cannot be empty.");

    auto localQ = settings.local();
    auto similarityRadius = settings.radius();
    auto minPts = settings.minimalClusterSize();

    ClusterLabels result;
    if(localQ) {
        cluster.permuteRelevantElectronsToFront(samples_);
        DensityBasedScan<double, Cluster, ClusteringMetric::distance<Spatial, Local>> dbscan(cluster);
        result = dbscan.findClusters(similarityRadius, minPts);
    } else {
        DensityBasedScan<double, Cluster, ClusteringMetric::distance<Spatial, Global>> dbscan(cluster);
        result = dbscan.findClusters(similarityRadius, minPts);
    }

    Cluster supercluster(static_cast<Cluster::size_type>(result.numberOfClusters));

    for (int i = 0; i < result.numberOfClusters; ++i)
        for (auto  [j, g] : enumerate(cluster))
            if (result.labels[j] == i)
                supercluster[i].emplace_back(std::move(g));

    orderByBestMatchDistance<Spatial>(supercluster, similarityRadius, localQ);

    cluster = supercluster;

    // sort by function value before leaving
    cluster.sortAll();
}

