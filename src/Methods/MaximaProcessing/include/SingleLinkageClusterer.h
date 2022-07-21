// Copyright (C) 2018-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SINGLELINKAGECLUSTERER_H
#define INPSIGHTS_SINGLELINKAGECLUSTERER_H

#include <DensityBasedScan.h>
#include <ISettings.h>
#include <spdlog/spdlog.h>
#include <IProcess.h>
#include <Cluster.h>
#include "ClusteringMetric.h"

namespace Settings {
    class SingleLinkageClusterer : public ISettings {
    public:
        Property<double> radius = {0.2, VARNAME(radius)};
        Property<size_t> minimalClusterSize = {1, VARNAME(minimalClusterSize)};
        Property<bool> local =  {false, VARNAME(local)}; // TODO unite all general clusterer settigns in a parent settings class
        Property<bool> sortRemainder = {false, VARNAME(sortRemainder)};

        SingleLinkageClusterer();
        explicit SingleLinkageClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::SingleLinkageClusterer)


class SingleLinkageClusterer : public IClusterer {
public:
    static Settings::SingleLinkageClusterer settings;

    explicit SingleLinkageClusterer(std::vector<Sample> &samples);

    void cluster(Cluster &group) override;

private:
    template<MetricType metric>
    void orderByBestMatchDistance(Cluster &supercluster, double threshold, bool localQ) const {
        for (auto &subcluster : supercluster) {
            sort(subcluster.begin(), subcluster.end());

            // starting sortedCluster with one active cluster, which is erased from subcluster
            Cluster sortedCluster({*subcluster.begin()});
            subcluster.erase(subcluster.begin());
            long activeClusters = 1;

            while (!subcluster.empty()) {
                // setting newClusters empty again
                Cluster newClusters;

                // iterating over all active clusters (at the end of sortedCluster)
                for (auto i = sortedCluster.end() - activeClusters; i != sortedCluster.end(); ++i) {
                    // iterating over all clusters remaining in the unsorted subcluster
                    for (auto j = subcluster.begin(); j != subcluster.end(); ++j) {
                        bool isSimilarQ = false;

                        if (localQ) {
                            if (i->getSelectedElectronsCount() == j->getSelectedElectronsCount())
                                isSimilarQ = compareAndPermuteLocal<metric>(threshold, *i, *j);
                        } else {
                            isSimilarQ = compareAndPermute<metric>(threshold, *i, *j);
                        }

                        if (isSimilarQ) {
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

    // TODO can these the permutation execution be shifted into the cluster via the bestMatchResult method?

    template<MetricType metric>
    bool compareAndPermute(double threshold, const Cluster &i, Cluster &j) const {
        bool isSimilarQ = false;

        auto[norm, perm] = ClusteringMetric::bestMatchResult<metric, Global>(j, i);

        if (norm <= threshold) {
            isSimilarQ = true;
            j.permuteAll(perm, samples_);
        };

        return isSimilarQ;
    }

    template<MetricType metric>
    bool compareAndPermuteLocal(double threshold, const Cluster &i, Cluster &j) const {

        auto electronsCount = i.representative()->maximum().numberOfEntities();
        auto iElectronsCount = i.getSelectedElectronsCount();
        auto jElectronsCount = j.getSelectedElectronsCount();

        auto[norm, perm] = ClusteringMetric::bestMatchResult<metric, Local>(j, i);

        bool isSimilarQ = false;

        if (norm <= threshold) {
            isSimilarQ = true;
            j.permuteAll(PermutationHandling::headToFullPermutation(perm, electronsCount), samples_);

            if (settings.sortRemainder()) { // TODO get rid of this
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
};

#endif //INPSIGHTS_SINGLELINKAGECLUSTERER_H
