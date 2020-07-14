/* Copyright (C) 2019 Leonard Reuter.
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

#include <ReferencePositionsClusterer.h>
#include <BestMatchDistance.h>
#include <BestMatch.h>
#include <ParticleSelection.h>
#include <Reference.h>
#include <Cluster.h>
#include <functional>
#include <ErrorHandling.h>

namespace Settings {
    ReferencePositionsClusterer::ReferencePositionsClusterer()
            : ISettings(VARNAME(ReferencePositionsClusterer)) {
        radius.onChange_.connect(
                [&](double value) {
                    if (value <= 0.0)
                        throw std::invalid_argument("The radius has to be larger than zero.");
                });

        local.onChange_.connect(
                [&](bool value) {
                    if (value) spdlog::info("Clustering is local");
                    else spdlog::info("Clustering is global.");
                });
    };

    ReferencePositionsClusterer::ReferencePositionsClusterer(const YAML::Node &node)
            : ReferencePositionsClusterer() {
        doubleProperty::decode(node, radius);
        boolProperty::decode(node, local);
        boolProperty::decode(node, sortRemainder);
    };

    void ReferencePositionsClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][local.name()] = local();
        node[className][sortRemainder.name()] = sortRemainder();
    };
}

YAML_SETTINGS_DEFINITION(Settings::ReferencePositionsClusterer)

Settings::ReferencePositionsClusterer ReferencePositionsClusterer::settings = Settings::ReferencePositionsClusterer();

ReferencePositionsClusterer::ReferencePositionsClusterer(std::vector<Sample> &samples)
        : IClusterer(samples)
          {}


void ReferencePositionsClusterer::cluster(Cluster &cluster) {
    assert(!cluster.empty() && "The cluster cannot be empty.");

    auto localQ = settings.local();
    auto similarityRadius = settings.radius();

    // sorting relevant electrons to the front
    if(localQ) cluster.permuteRelevantElectronsToFront(samples_);

    // initialize supercluster
    Cluster superCluster({Cluster({*cluster.begin()})});

    for (auto subCluster = std::next(cluster.begin()); subCluster != cluster.end(); ++subCluster) {
        // bool to decide, whether subCluster of cluster is added to a sortedCluster of superCluster
        // or to superCluster as a new sortedCluster
        bool isSimilarQ = false;
        for (auto sortedCluster = superCluster.begin(); sortedCluster != superCluster.end(); ++sortedCluster) {

            if(localQ){
                if(subCluster->getSelectedElectronsCount() == sortedCluster->getSelectedElectronsCount())
                    isSimilarQ = compareLocal(sortedCluster, subCluster, similarityRadius);
            } else {
                isSimilarQ = compareGlobal(sortedCluster, subCluster, similarityRadius);
            }

            if (isSimilarQ) break;
        }
        // if loop ended and no similar cluster was found
        if (!isSimilarQ) {
            // adds subCluster as a new sortedCluster to superCluster
            superCluster.emplace_back(Cluster({*subCluster}));
        }
    }
    cluster = superCluster;

    // sort by function value before leaving
    cluster.sort();
}

bool ReferencePositionsClusterer::compareLocal(std::vector<Cluster>::iterator &sortedCluster, std::vector<Cluster>::iterator &subCluster,
                                     double similarityRadius) const {
    bool isSimilarQ = false;

    // only check similarity of sortedCluster and subCluster, if the number of selected indices is equal
    // this requires sortedCluster having the correct electrons count
    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
            subCluster->representative()->maximum().head(
                    subCluster->getSelectedElectronsCount()).positionsVector(),
            sortedCluster->representative()->maximum().head(
                    sortedCluster->getSelectedElectronsCount()).positionsVector());

    if (norm < similarityRadius) {
        isSimilarQ = true;

        auto  electronsNumber = sortedCluster->representative()->maximum().numberOfEntities();
        subCluster->permuteAll(BestMatch::headToFullPermutation(perm, electronsNumber), samples_);
        if (settings.sortRemainder()) { //TODO remove dependency
            auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subCluster->representative()->maximum().tail(
                            electronsNumber - subCluster->getSelectedElectronsCount()).positionsVector(),
                    sortedCluster->representative()->maximum().tail(
                            electronsNumber - sortedCluster->getSelectedElectronsCount()).positionsVector());
            subCluster->permuteAll(BestMatch::tailToFullPermutation(perm, electronsNumber), samples_);
        }
        sortedCluster->emplace_back(*subCluster);
    }

    return isSimilarQ;
}

bool ReferencePositionsClusterer::compareGlobal(std::vector<Cluster>::iterator &sortedCluster, std::vector<Cluster>::iterator &subCluster,
                                               double similarityRadius) const {
    bool isSimilarQ = false;

    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
            subCluster->representative()->maximum().positionsVector(),
            sortedCluster->representative()->maximum().positionsVector());

    if (norm < similarityRadius) {
        isSimilarQ = true;

        subCluster->permuteAll(perm, samples_);
        sortedCluster->emplace_back(*subCluster);
    }

    return isSimilarQ;
}
