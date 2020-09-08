// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <IdentityClusterer.h>
#include <PreClusterer.h>
#include <Reference.h>
#include <BestMatchDistance.h>

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
        
// assumes a sorted reference vector
void PreClusterer::cluster(Cluster& cluster) {
    assert(!cluster.empty() && "The cluster cannot be empty.");

    auto similarityRadius = settings.radius();
    auto valueIncrement = settings.valueIncrement();
    auto atoms = cluster.representative()->nuclei();

    // insert first element
    Cluster supercluster({Cluster({*cluster.begin()})});
    cluster.erase(cluster.begin());

    //Presort
    for (auto subcluster = cluster.begin(); subcluster != cluster.end(); ++subcluster) {

        Cluster lowerRef(Reference(atoms, subcluster->representative()->value() - valueIncrement,
                ElectronsVector(), 0));
        Cluster upperRef(Reference(atoms, subcluster->representative()->value() + valueIncrement,
                ElectronsVector(), 0));

        auto superclusterLowerBoundIt = std::lower_bound(
                supercluster.begin(),
                supercluster.end(),
                lowerRef);
        auto superclusterUpperBoundIt = std::upper_bound(
                supercluster.begin(),
                supercluster.end(),
                upperRef);

        // iterate over all supercluster members within the value range
        std::list<bool> outsideQ;
        for(auto subclusterFromSuperclusterBoundaries = superclusterLowerBoundIt;
            subclusterFromSuperclusterBoundaries != superclusterUpperBoundIt; ++subclusterFromSuperclusterBoundaries) {

            auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subcluster->representative()->maximum().positionsVector(),
                    subclusterFromSuperclusterBoundaries->representative()->maximum().positionsVector());
            if(norm > similarityRadius)
                outsideQ.emplace_back(true);
            else
                outsideQ.emplace_back(false);
        }
        if(std::all_of(outsideQ.begin(), outsideQ.end(), [](bool b){return b;})) {
            supercluster.emplace_back(Cluster({*subcluster}));
            subcluster = cluster.erase(subcluster);
            --subcluster;
        }
    }

    // start with the second subcluster
    for (auto subcluster = cluster.begin(); subcluster != cluster.end(); ++subcluster) {

        // Define value range of the supercluster
        Cluster lowerRef(Reference(atoms, subcluster->representative()->value() - valueIncrement,
                ElectronsVector(), 0));
        Cluster upperRef(Reference(atoms, subcluster->representative()->value() + valueIncrement,
                ElectronsVector(), 0));

        auto superclusterLowerBoundIt = std::lower_bound(
                supercluster.begin(),
                supercluster.end(),
                lowerRef);
        auto superclusterUpperBoundIt = std::upper_bound(
                supercluster.begin(),
                supercluster.end(),
                upperRef);

        // iterate over all supercluster members within the value range
        auto overallBestMatchNorm = std::numeric_limits<double>::max();
        auto overallBestMatchPerm =
                Eigen::PermutationMatrix<Eigen::Dynamic>(cluster.representative()->maximum().numberOfEntities());
        auto bestMatchSubclusterFromSuperclusterBoundaries = superclusterLowerBoundIt;
        for (auto subclusterFromSuperclusterBoundaries = superclusterLowerBoundIt;
        subclusterFromSuperclusterBoundaries != superclusterUpperBoundIt; ++subclusterFromSuperclusterBoundaries) {

            auto [norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subcluster->representative()->maximum().positionsVector(),
                    subclusterFromSuperclusterBoundaries->representative()->maximum().positionsVector());

            if (norm <= overallBestMatchNorm) {
                overallBestMatchNorm = norm;
                overallBestMatchPerm = perm;
                bestMatchSubclusterFromSuperclusterBoundaries = subclusterFromSuperclusterBoundaries;
            }

        }
        if (overallBestMatchNorm <= similarityRadius) {
            subcluster->permuteAll(overallBestMatchPerm, samples_);
            bestMatchSubclusterFromSuperclusterBoundaries->emplace_back(*subcluster);
        }
        else {
            throw std::exception();
        }
    }
    cluster = supercluster;

    // sort by function value before leaving
    cluster.sortAll();
}
