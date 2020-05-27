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
#include <Group.h>
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


void ReferencePositionsClusterer::cluster(Group &group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto localQ = settings.local();
    auto similarityRadius = settings.radius();

    // sorting relevant electrons to the front
    if(localQ) group.permuteRelevantElectronsToFront(samples_);

    // initialize supergroup
    Group superGroup({Group({*group.begin()})});

    for (auto subGroup = std::next(group.begin()); subGroup != group.end(); ++subGroup) {
        // bool to decide, whether subGroup of group is added to a sortedGroup of superGroup
        // or to superGroup as a new sortedGroup
        bool isSimilarQ = false;
        for (auto sortedGroup = superGroup.begin(); sortedGroup != superGroup.end(); ++sortedGroup) {

            if(localQ){
                if(subGroup->getSelectedElectronsCount() == sortedGroup->getSelectedElectronsCount())
                    isSimilarQ = compareLocal(sortedGroup, subGroup, similarityRadius);
            } else {
                isSimilarQ = compareGlobal(sortedGroup, subGroup, similarityRadius);
            }

            if (isSimilarQ) break;
        }
        // if loop ended and no similar group was found
        if (!isSimilarQ) {
            // adds subGroup as a new sortedGroup to superGroup
            superGroup.emplace_back(Group({*subGroup}));
        }
    }
    group = superGroup;

    // sort by function value before leaving
    group.sort();
}

bool ReferencePositionsClusterer::compareLocal(std::vector<Group>::iterator &sortedGroup, std::vector<Group>::iterator &subGroup,
                                     double similarityRadius) const {
    bool isSimilarQ = false;

    // only check similarity of sortedGroup and subGroup, if the number of selected indices is equal
    // this requires sortedGroup having the correct electrons count
    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
            subGroup->representative()->maximum().head(
                    subGroup->getSelectedElectronsCount()).positionsVector(),
            sortedGroup->representative()->maximum().head(
                    sortedGroup->getSelectedElectronsCount()).positionsVector());

    if (norm < similarityRadius) {
        isSimilarQ = true;

        auto  electronsNumber = sortedGroup->representative()->maximum().numberOfEntities();
        subGroup->permuteAll(BestMatch::headToFullPermutation(perm, electronsNumber), samples_);
        if (settings.sortRemainder()) { //TODO remove dependency
            auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subGroup->representative()->maximum().tail(
                            electronsNumber - subGroup->getSelectedElectronsCount()).positionsVector(),
                    sortedGroup->representative()->maximum().tail(
                            electronsNumber - sortedGroup->getSelectedElectronsCount()).positionsVector());
            subGroup->permuteAll(BestMatch::tailToFullPermutation(perm, electronsNumber), samples_);
        }
        sortedGroup->emplace_back(*subGroup);
    }

    return isSimilarQ;
}

bool ReferencePositionsClusterer::compareGlobal(std::vector<Group>::iterator &sortedGroup, std::vector<Group>::iterator &subGroup,
                                               double similarityRadius) const {
    bool isSimilarQ = false;

    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
            subGroup->representative()->maximum().positionsVector(),
            sortedGroup->representative()->maximum().positionsVector());

    if (norm < similarityRadius) {
        isSimilarQ = true;

        subGroup->permuteAll(perm, samples_);
        sortedGroup->emplace_back(*subGroup);
    }

    return isSimilarQ;
}
