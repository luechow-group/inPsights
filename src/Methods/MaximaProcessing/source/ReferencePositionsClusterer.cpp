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
#include <NearestElectrons.h>
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
    };

    ReferencePositionsClusterer::ReferencePositionsClusterer(const YAML::Node &node)
            : ReferencePositionsClusterer() {
        doubleProperty::decode(node, radius);
    };

    void ReferencePositionsClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
    };
}

YAML_SETTINGS_DEFINITION(Settings::ReferencePositionsClusterer)

Settings::ReferencePositionsClusterer ReferencePositionsClusterer::settings = Settings::ReferencePositionsClusterer();

ReferencePositionsClusterer::ReferencePositionsClusterer(std::vector<Sample> &samples)
        : IClusterer(samples)
          {}

void ReferencePositionsClusterer::cluster(Group &group) {
    assert(!group.empty() && "The group cannot be empty.");

    // sorting relevant electrons to the front
    group.permuteRelevantElectronsToFront(samples_);

    // initialize supergroup
    Group superGroup({Group({*group.begin()})});

    // bool to decide, whether subGroup of group is added to a sortedGroup of superGroup
    // or to superGroup as a new sortedGroup
    bool isSimilarQ;

    auto similarityRadius = settings.radius();
    long electronsNumber = group.representative()->maximum().numberOfEntities();

    for (auto subGroup = std::next(group.begin()); subGroup != group.end(); ++subGroup) {
        isSimilarQ = false;
        for (auto sortedGroup = superGroup.begin(); sortedGroup != superGroup.end(); ++sortedGroup) {
            // only check similarity of sortedGroup and subGroup, if the number of selected indices is equal
            if (subGroup->getSelectedElectronsCount() == superGroup.getSelectedElectronsCount()) {
                auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                        subGroup->representative()->maximum().getFirstParticles(
                                subGroup->getSelectedElectronsCount()).positionsVector(),
                        sortedGroup->representative()->maximum().getFirstParticles(
                                superGroup.getSelectedElectronsCount()).positionsVector());

                if (norm < similarityRadius) {
                    subGroup->permuteAll(BestMatch::headToFullPermutation(perm, electronsNumber), samples_);
                    if (settings.sortRemainder()) {
                        auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                                subGroup->representative()->maximum().tail(
                                        electronsNumber - subGroup->getSelectedElectronsCount()).positionsVector(),
                                sortedGroup->representative()->maximum().tail(
                                        electronsNumber - superGroup.getSelectedElectronsCount()).positionsVector());
                        subGroup->permuteAll(BestMatch::tailToFullPermutation(perm, electronsNumber), samples_);
                    }
                    sortedGroup->emplace_back(*subGroup);
                    isSimilarQ = true;
                    break;
                }
            }
        }
        if (!isSimilarQ) {
            // adds subGroup as a new sortedGroup to superGroup
            superGroup.emplace_back(Group({*subGroup}));
        }
    }
    group = superGroup;

    // sort by function value before leaving
    group.sort();
}
