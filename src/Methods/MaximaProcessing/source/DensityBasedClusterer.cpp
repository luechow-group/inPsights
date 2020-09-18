// Copyright (C) 2018-2020 Michael Heuer.
// Copyright (C) 2018 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <DensityBasedClusterer.h>
#include <PreClusterer.h>
#include <DistanceBasedMetric.h>
#include <Enumerate.h>
#include <Reference.h>
#include "ClusteringMetric.h"

namespace Settings {
    DensityBasedClusterer::DensityBasedClusterer()
    : ISettings(VARNAME(DensityBasedClusterer)) {
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

    DensityBasedClusterer::DensityBasedClusterer(const YAML::Node &node)
            : DensityBasedClusterer() {
        doubleProperty::decode(node, radius);
        size_tProperty ::decode(node, minimalClusterSize);
        boolProperty::decode(node, local);
        boolProperty::decode(node, sortRemainder);
    }

    void DensityBasedClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][minimalClusterSize.name()] = minimalClusterSize();
        node[className][local.name()] = local();
        node[className][sortRemainder.name()] = sortRemainder();
    }
}
YAML_SETTINGS_DEFINITION(Settings::DensityBasedClusterer)

Settings::DensityBasedClusterer DensityBasedClusterer::settings = Settings::DensityBasedClusterer();


DensityBasedClusterer::DensityBasedClusterer(std::vector<Sample> &samples)
        : IClusterer(samples) {};

void DensityBasedClusterer::cluster(Cluster& cluster) {
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

    orderByBestMatchDistance(supercluster, similarityRadius, localQ);

    cluster = supercluster;

    // sort by function value before leaving
    cluster.sortAll();
}

void DensityBasedClusterer::orderByBestMatchDistance(Cluster &supercluster, double threshold, bool localQ) const {
    for (auto &subcluster : supercluster) {
        sort(subcluster.begin(), subcluster.end());

        // starting sortedCluster with one active cluster, which is erased from subcluster
        Cluster sortedCluster({*subcluster.begin()});
        subcluster.erase(subcluster.begin());
        long activeClusters = 1;

        while (!subcluster.empty()){
            // setting newClusters empty again
            Cluster newClusters;

            // iterating over all active clusters (at the end of sortedCluster)
            for (auto i = sortedCluster.end() - activeClusters; i != sortedCluster.end(); ++i){
                // iterating over all clusters remaining in the unsorted subcluster
                for (auto j = subcluster.begin(); j != subcluster.end(); ++j) {
                    bool isSimilarQ = false;

                    if(localQ) {
                        if (i->getSelectedElectronsCount() == j->getSelectedElectronsCount())
                            isSimilarQ = compareLocal(threshold, *i, *j);
                    } else {
                        isSimilarQ = compareAndPermute(threshold, *i, *j);
                    }
                    
                    if(isSimilarQ) {
                        // moving j from subcluster to newClusters
                        newClusters.emplace_back(*j);
                        j = subcluster.erase(j);

                        // the iterator has to be set back by one because the j element was erased and
                        // ++j of the for loop would otherwise skip one cluster of subcluster
                        --j;
                    }
                }
            }
            // moving all clusters from newClusters to sortedCluster()
            activeClusters = newClusters.size();
            for (auto &newCluster : newClusters) {
                sortedCluster.emplace_back(newCluster);
            };
        };

        subcluster = sortedCluster;
    }
}

bool DensityBasedClusterer::compareAndPermute(double threshold, const Cluster &i, Cluster &j) const {
    bool isSimilarQ = false;

    auto[norm, perm] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
                    j.representative()->maximum().positionsVector(),
                    i.representative()->maximum().positionsVector());

    if (norm <= threshold) {
        isSimilarQ = true;
        j.permuteAll(perm, samples_);
    };

    return isSimilarQ;
}

bool DensityBasedClusterer::compareLocal(double threshold, const Cluster &i, Cluster &j) const {

    auto electronsCount = i.representative()->maximum().numberOfEntities();
    auto iElectronsCount = i.getSelectedElectronsCount();
    auto jElectronsCount = j.getSelectedElectronsCount();

    auto[norm, perm] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
            j.representative()->maximum().head(jElectronsCount).positionsVector(),
            i.representative()->maximum().head(iElectronsCount).positionsVector());

    bool isSimilarQ = false;

    if (norm <= threshold) {
        isSimilarQ = true;
        j.permuteAll(PermutationHandling::headToFullPermutation(perm, electronsCount), samples_);
        if (settings.sortRemainder()) {
            auto[norm, perm] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
                    j.representative()->maximum().tail(
                            electronsCount - jElectronsCount).positionsVector(),
                    i.representative()->maximum().tail(
                            electronsCount - iElectronsCount).positionsVector());
            j.permuteAll(PermutationHandling::tailToFullPermutation(perm, electronsCount),
                          samples_);
        }
    };

    return isSimilarQ;
}
