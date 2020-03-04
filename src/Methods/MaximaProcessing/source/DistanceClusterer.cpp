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
#include <DistanceClusterer.h>
#include <BestMatchDistance.h>
#include <ValueSorter.h>

namespace Settings {
    DistanceClusterer::DistanceClusterer()
    : ISettings(VARNAME(DistanceClusterer)) {}

    DistanceClusterer::DistanceClusterer(const YAML::Node &node)
            : DistanceClusterer() {
        doubleProperty::decode(node, radius);
        doubleProperty::decode(node, valueIncrement);
    }

    void DistanceClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][valueIncrement.name()] = valueIncrement();
    }
}
YAML_SETTINGS_DEFINITION(Settings::DistanceClusterer)

Settings::DistanceClusterer DistanceClusterer::settings = Settings::DistanceClusterer();


DistanceClusterer::DistanceClusterer(std::vector<Sample> &samples)
        : IClusterer(samples){}
        
// assumes a sorted reference vector
void DistanceClusterer::cluster(Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto similarityRadius = settings.radius();
    auto valueIncrement = settings.valueIncrement();

    // insert first element
    Group supergroup({Group({*group.begin()})});
    group.erase(group.begin());

    //Presort
    for (auto subgroup = group.begin(); subgroup != group.end(); ++subgroup) {

        Group lowerRef(Reference(subgroup->representative()->value() - valueIncrement));
        Group upperRef(Reference(subgroup->representative()->value() + valueIncrement));

        auto supergroupLowerBoundIt = std::lower_bound(
                supergroup.begin(),
                supergroup.end(),
                lowerRef);
        auto supergroupUpperBoundIt = std::upper_bound(
                supergroup.begin(),
                supergroup.end(),
                upperRef);

        // iterate over all supergroup members within the value range
        std::list<bool> outsideQ;
        for(auto subgroupFromSupergroupBoundaries = supergroupLowerBoundIt;
            subgroupFromSupergroupBoundaries != supergroupUpperBoundIt; ++subgroupFromSupergroupBoundaries) {

            auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subgroup->representative()->maximum().positionsVector(),
                    subgroupFromSupergroupBoundaries->representative()->maximum().positionsVector());
            if(norm > similarityRadius)
                outsideQ.emplace_back(true);
            else
                outsideQ.emplace_back(false);
        }
        if(std::all_of(outsideQ.begin(), outsideQ.end(), [](bool b){return b;})) {
            supergroup.emplace_back(Group({*subgroup}));
            subgroup = group.erase(subgroup);
            --subgroup;
        }
    }

    // start with the second subgroup
    for (auto subgroup = group.begin(); subgroup != group.end(); ++subgroup) {

        // Define value range of the supergroup
        Group lowerRef(Reference(subgroup->representative()->value() - valueIncrement));
        Group upperRef(Reference(subgroup->representative()->value() + valueIncrement));

        auto supergroupLowerBoundIt = std::lower_bound(
                supergroup.begin(),
                supergroup.end(),
                lowerRef);
        auto supergroupUpperBoundIt = std::upper_bound(
                supergroup.begin(),
                supergroup.end(),
                upperRef);

        // iterate over all supergroup members within the value range
        auto overallBestMatchNorm = std::numeric_limits<double>::max();
        auto overallBestMatchPerm =
                Eigen::PermutationMatrix<Eigen::Dynamic>(group.representative()->maximum().numberOfEntities());
        auto bestMatchSubgroupFromSupergroupBoundaries = supergroupLowerBoundIt;
        for (auto subgroupFromSupergroupBoundaries = supergroupLowerBoundIt;
        subgroupFromSupergroupBoundaries != supergroupUpperBoundIt; ++subgroupFromSupergroupBoundaries) {

            auto [norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subgroup->representative()->maximum().positionsVector(),
                    subgroupFromSupergroupBoundaries->representative()->maximum().positionsVector());

            if (norm <= overallBestMatchNorm) {
                overallBestMatchNorm = norm;
                overallBestMatchPerm = perm;
                bestMatchSubgroupFromSupergroupBoundaries = subgroupFromSupergroupBoundaries;
            }

        }
        if (overallBestMatchNorm <= similarityRadius) {
            subgroup->permuteAll(overallBestMatchPerm, samples_);
            bestMatchSubgroupFromSupergroupBoundaries->emplace_back(*subgroup);
        }
        else {
            throw std::exception();
        }
    }
    group = supergroup;

    // sort by function value before leaving
    group.sortAll();
}
