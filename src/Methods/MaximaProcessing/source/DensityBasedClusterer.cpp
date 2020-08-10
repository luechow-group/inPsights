/* Copyright (C) 2018-2019 Michael Heuer.
 * Copyright (C) 2018 Leonard Reuter.
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

#include <DensityBasedClusterer.h>
#include <PreClusterer.h>
#include <BestMatchDistance.h>
#include <Enumerate.h>
#include <Reference.h>

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

double DensityBasedClusterer::wrapper(const Cluster &g1, const Cluster &g2) {
    return BestMatch::Distance::compare<Eigen::Infinity, 2>(
            g1.representative()->maximum().positionsVector(),
            g2.representative()->maximum().positionsVector()).metric;
};

double DensityBasedClusterer::wrapperLocal(const Cluster &g1, const Cluster &g2) {

    auto g1ElectronsCount = g1.getSelectedElectronsCount();
    auto g2ElectronsCount = g2.getSelectedElectronsCount();

    if (g1ElectronsCount == g2ElectronsCount) {
        return BestMatch::Distance::compare<Eigen::Infinity, 2>(
                g1.representative()->maximum().head(g1ElectronsCount).positionsVector(),
                g2.representative()->maximum().head(g2ElectronsCount).positionsVector()).metric;
    }

    return std::numeric_limits<double>::max();
};

void DensityBasedClusterer::cluster(Cluster& cluster) {
    assert(!cluster.empty() && "The cluster cannot be empty.");

    auto localQ = settings.local();
    auto similarityRadius = settings.radius();
    auto minPts = settings.minimalClusterSize();

    ClusterLabels result;
    if(localQ) {
        cluster.permuteRelevantElectronsToFront(samples_);
        DensityBasedScan<double, Cluster, DensityBasedClusterer::wrapperLocal> dbscan(cluster);
        result = dbscan.findClusters(similarityRadius, minPts);
    } else {
        DensityBasedScan<double, Cluster, DensityBasedClusterer::wrapper> dbscan(cluster);
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
                            isSimilarQ = compareLocal(threshold, subcluster, newClusters, i, j);
                    } else {
                        isSimilarQ = compare(threshold, subcluster, newClusters, i, j);
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

bool DensityBasedClusterer::compare(double threshold, Cluster &subcluster, Cluster &newClusters,
        const std::vector<Cluster>::iterator &i,
        std::vector<Cluster>::iterator &j) const {
    bool isSimilarQ = false;

    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    j->representative()->maximum().positionsVector(),
                    i->representative()->maximum().positionsVector());

    if (norm <= threshold) {
        isSimilarQ = true;
        j->permuteAll(perm, samples_);
    };

    return isSimilarQ;
}

bool DensityBasedClusterer::compareLocal(double threshold, Cluster &subcluster, Cluster &newClusters,
                                         const std::vector<Cluster>::iterator &i,
                                         std::vector<Cluster>::iterator &j) const {
    auto electronsCount = subcluster.representative()->maximum().numberOfEntities();
    auto iElectronsCount = i->getSelectedElectronsCount();
    auto jElectronsCount = j->getSelectedElectronsCount();

    bool isSimilarQ = false;

    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
            j->representative()->maximum().head(jElectronsCount).positionsVector(),
            i->representative()->maximum().head(iElectronsCount).positionsVector());

    if (norm <= threshold) {
        isSimilarQ = true;
        j->permuteAll(BestMatch::headToFullPermutation(perm, electronsCount), samples_);
        if (settings.sortRemainder()) {
            auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    j->representative()->maximum().tail(
                            electronsCount - jElectronsCount).positionsVector(),
                    i->representative()->maximum().tail(
                            electronsCount - iElectronsCount).positionsVector());
            j->permuteAll(BestMatch::tailToFullPermutation(perm, electronsCount),
                          samples_);
        }
    };

    return isSimilarQ;
}
