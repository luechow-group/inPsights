// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <IdentityClusterer.h>
#include <PreClusterer.h>
#include <Maximum.h>
#include <DistanceBasedMetric.h>

namespace Settings {
    PreClusterer::PreClusterer()
    : ISettings(VARNAME(PreClusterer)) {}

    PreClusterer::PreClusterer(const YAML::Node &node)
            : PreClusterer() {
        doubleProperty::decode(node, radius);
        doubleProperty::decode(node, valueIncrement);
    }

    void PreClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][valueIncrement.name()] = valueIncrement();
    }
}
YAML_SETTINGS_DEFINITION(Settings::PreClusterer)

Settings::PreClusterer PreClusterer::settings = Settings::PreClusterer();


PreClusterer::PreClusterer(std::vector<Sample> &samples)
        : IClusterer(samples){}
        
// assumes a sorted maxima vector
void PreClusterer::cluster(Cluster& cluster) {
    assert(!cluster.empty() && "The cluster cannot be empty.");

    auto similarityRadius = settings.radius();
    auto valueIncrement = settings.valueIncrement();
    auto atoms = cluster.representative()->nuclei();

    // The first cluster becomes the initial supercluster
    Cluster superclusters({Cluster({*cluster.begin()})});
    cluster.erase(cluster.begin());

    // Iterate over all remaining subclusters (this update on erase)
    for (auto subcluster = cluster.begin(); subcluster != cluster.end(); ++subcluster) {

        // Define function value window
        Cluster lowerRef(Maximum(atoms, subcluster->representative()->value() - valueIncrement,
                ElectronsVector(), 0));
        Cluster upperRef(Maximum(atoms, subcluster->representative()->value() + valueIncrement,
                ElectronsVector(), 0));

        auto superclustersLowerBoundIt = std::lower_bound(
                superclusters.begin(),
                superclusters.end(),
                lowerRef);
        auto superclustersUpperBoundIt = std::upper_bound(
                superclusters.begin(),
                superclusters.end(),
                upperRef);

        // iterate over all supercluster members within the value range and check,
        // if the current subcluster falls inside the similarity radius of any supercluster
        bool insideAnySuperclustersQ = false;
        for(auto superclusterFromBoundaries = superclustersLowerBoundIt;
            superclusterFromBoundaries != superclustersUpperBoundIt; ++superclusterFromBoundaries) {

            // comparison
            auto[norm, perm] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
                    subcluster->representative()->maximum().positionsVector(),
                    superclusterFromBoundaries->representative()->maximum().positionsVector());

            // if inside the similarity radius of the current supercluster, put it in there
            // else, make notice in outsideQ
            if(norm <= similarityRadius) {
                subcluster->permuteAll(perm, samples_);
                superclusterFromBoundaries->emplace_back(*subcluster);
                insideAnySuperclustersQ = true;
                break;
            }
        }

        // if the subcluster does not fall inside any supercluster, make it a supercluster itself
        if(!insideAnySuperclustersQ) {

            // emplace a copy in the supercluster
            superclusters.emplace_back(Cluster({*subcluster}));

            // erase the subcluster from the list of subclusters
            subcluster = cluster.erase(subcluster);
            --subcluster;
        }
    }

    cluster = superclusters;

    // sort by function value before leaving
    cluster.sortAll();
}
