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

#include <IdentityClusterer.h>
#include <PreClusterer.h>
#include <BestMatchDistance.h>
#include <ValueSorter.h>
#include <spdlog/spdlog.h>

namespace Settings {
    IdentityClusterer::IdentityClusterer()
    : ISettings(VARNAME(IdentityClusterer)) {}

    IdentityClusterer::IdentityClusterer(const YAML::Node &node)
            : IdentityClusterer() {
        doubleProperty::decode(node, radius);
        doubleProperty::decode(node, valueIncrement);
    }

    void IdentityClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][valueIncrement.name()] = valueIncrement();
    }
}
YAML_SETTINGS_DEFINITION(Settings::IdentityClusterer)

Settings::IdentityClusterer IdentityClusterer::settings = Settings::IdentityClusterer();


IdentityClusterer::IdentityClusterer(std::vector<Sample> &samples)
        : IClusterer(samples) {}

void IdentityClusterer::cluster(Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto identityRadius = settings.radius();
    auto valueIncrement = settings.valueIncrement();
    auto atoms = group.representative()->nuclei();

    auto beginIt = group.begin();

    while (beginIt != group.end()) {
        auto total = group.size();//std::distance(group.begin(), group.end());
        auto endIt = std::upper_bound(beginIt, group.end(), Group(Reference(atoms,
                                                                            beginIt->representative()->value() +
                                                                            valueIncrement, ElectronsVector(), 0)));

        spdlog::info("Global identiy search in interval {} to {}, total: {}",
                      total - std::distance(beginIt, group.end()),
                      total - std::distance(endIt, group.end()),
                      std::distance(group.begin(), group.end()));

        auto it = beginIt;

        if (beginIt != endIt) {
            it++; // start with the element next to beginIt
            while (it != endIt)
                subLoop(group, beginIt, it, endIt, identityRadius, valueIncrement);

            beginIt = endIt;
        } else ++beginIt; // range is zero
    }
    // sort by function value before leaving
    group.sortAll();
}

void IdentityClusterer::subLoop(Group& group,
        Group::iterator &beginIt,
        Group::iterator &it,
        Group::iterator &endIt,
        double distThresh,
        double valueIncrement) {

    auto atoms = group.representative()->nuclei();

    //TODO calculate only alpha electron distances and skip beta electron hungarian if dist is too large
    auto [norm, perm] = BestMatch::Distance::compare<Spin, Eigen::Infinity, 2>(
            it->representative()->maximum(),
            (*beginIt).representative()->maximum());

    if (beginIt->representative()->maximum().typesVector().multiplicity() == 1) { // consider spin flip

        auto permuteeSpinFlipped = it->representative()->maximum();
        permuteeSpinFlipped.typesVector().flipSpins();

        auto [normFlipped, permFlipped] =
        BestMatch::Distance::compare<Spin, Eigen::Infinity, 2>(permuteeSpinFlipped, beginIt->representative()->maximum());

        if ((norm <= distThresh) || (normFlipped <= distThresh)) {
            if (norm <= normFlipped)
                addReference(group, beginIt, it, perm);
            else
                addReference(group, beginIt, it, permFlipped);
            endIt = std::upper_bound(beginIt, group.end(),
                    Group(Reference(atoms, beginIt->representative()->value() + valueIncrement,
                                    ElectronsVector(), 0)));
        } else it++;
    } else {  // don't consider spin flip
        if (norm <= distThresh) {
            addReference(group, beginIt, it, perm);
            endIt = std::upper_bound(beginIt, group.end(),
                    Group(Reference(atoms, beginIt->representative()->value() + valueIncrement,
                                    ElectronsVector(), 0)));
        } else it++;
    }
}

// TODO This method should be located inside of a reference container class
void IdentityClusterer::addReference(Group& group,
        const Group::iterator &beginIt,
        Group::iterator &it,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const {

    it->permuteAll(bestMatch, samples_);
    beginIt->representative()->mergeReference(it);
    it = group.erase(it); // erase returns the iterator of the following element
}