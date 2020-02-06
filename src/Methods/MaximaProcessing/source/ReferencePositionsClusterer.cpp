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

namespace Settings {
    ReferencePositionsClusterer::ReferencePositionsClusterer()
            : ISettings(VARNAME(ReferencePositionsClusterer)) {
        distanceMode.onChange_.connect(
                [&](std::string value) {
                    if (not (value == "minimum" || value == "average"))
                        throw std::invalid_argument("The distanceMode has to be minimum or average.");
                });
        radius.onChange_.connect(
                [&](double value) {
                    if (not (value > 0.0))
                        throw std::invalid_argument("The radius has to be larger than zero.");
                });
        maximalDistance.onChange_.connect(
                [&](double value) {
                    if (not (value > 0.0))
                        throw std::invalid_argument("The maximalDistance has to be larger than zero.");
                });
        maximalCount.onChange_.connect(
                [&](long value) {
                    if (not (value > 0))
                        throw std::invalid_argument("The maximalCount has to be larger than zero.");
                });
    };

    ReferencePositionsClusterer::ReferencePositionsClusterer(const YAML::Node &node)
            : ReferencePositionsClusterer() {
        doubleProperty::decode(node, radius);
        doubleProperty::decode(node, maximalDistance);
        longProperty::decode(node, maximalCount);
        stringProperty::decode(node, distanceMode);
        boolProperty::decode(node, invertSelection);
        boolProperty::decode(node, valenceOnly);
        boolProperty::decode(node, sortRemainder);
    };

    void ReferencePositionsClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][maximalDistance.name()] = maximalDistance();
        node[className][maximalCount.name()] = maximalCount();
        node[className][distanceMode.name()] = distanceMode();
        node[className][invertSelection.name()] = invertSelection();
        node[className][valenceOnly.name()] = valenceOnly();
        node[className][sortRemainder.name()] = sortRemainder();
    };
}

YAML_SETTINGS_DEFINITION(Settings::ReferencePositionsClusterer)

Settings::ReferencePositionsClusterer ReferencePositionsClusterer::settings = Settings::ReferencePositionsClusterer();

ReferencePositionsClusterer::ReferencePositionsClusterer(std::vector<Sample> &samples, AtomsVector &nuclei,
        std::vector<Eigen::Vector3d> &positions)
        : samples_(samples),
          nuclei_(nuclei),
          positions_(positions){
    if (settings.distanceMode() == "average"){
        distanceFunction_ = Metrics::averageDistance<2>;
    }
    else if (settings.distanceMode() == "minimum"){
        distanceFunction_ = Metrics::minimalDistance<2>;
    }
}

void ReferencePositionsClusterer::cluster(Group &group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto similarityRadius = settings.radius();
    long electronsNumber = group.representative()->maximum().numberOfEntities();

    std::list<long> subIndices;
    Eigen::PermutationMatrix<Eigen::Dynamic> permutation;

    // for every 'subGroup' in 'group', the number of relevant electrons is stored in 'counts'
    std::vector<long> counts;

    // sorting relevant electrons to the front
    for (auto subGroup = group.begin(); subGroup != group.end(); ++subGroup) {
        subIndices = ReferencePositionsClusterer::getRelevantIndices(subGroup->representative()->maximum());
        if (not settings.invertSelection()){
            counts.emplace_back(subIndices.size());
            // sorting will take the front indices, so they have to be permuted to the front
            permutation = BestMatch::getPermutationToFront(subIndices, electronsNumber);
        }
        else{
            if (settings.valenceOnly()){
                // since selection should be inverted, core indices have to be added to 'subIndices' before inverting
                subIndices.splice(subIndices.end(),
                        NearestElectrons::getNonValenceIndices(subGroup->representative()->maximum(), nuclei_));
            }
            counts.emplace_back(electronsNumber - subIndices.size());
            // sorting will take the front indices. Since the invertSelection is true,
            // 'subIndices' are permuted to the back
            permutation = BestMatch::getPermutationToBack(subIndices, electronsNumber);
        }
        subGroup->permuteAll(permutation, samples_);
    }

    std::vector<long> countsSuperGroup;
    auto countIterator = counts.begin();

    Group superGroup({Group({*group.begin()})});
    countsSuperGroup.emplace_back(*countIterator);
    countIterator++;

    auto countIteratorSuperGroup = countsSuperGroup.begin();

    // bool to decide, whether subGroup of group is added to a sortedGroup of superGroup
    // or to superGroup as a new sortedGroup
    bool isSimilarQ;

    for (auto subGroup = std::next(group.begin()); subGroup != group.end(); ++subGroup) {
        isSimilarQ = false;
        for (auto sortedGroup = superGroup.begin(); sortedGroup != superGroup.end(); ++sortedGroup) {
            // only check similarity of sortedGroup and subGroup, if the number of relevant indices ('count') is equal
            if (*countIterator == *countIteratorSuperGroup){
                auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                        subGroup->representative()->maximum().getFirstParticles(*countIterator).positionsVector(),
                        sortedGroup->representative()->maximum().getFirstParticles(*countIteratorSuperGroup).positionsVector());

                if (norm < similarityRadius) {
                    subGroup->permuteAll(BestMatch::headToFullPermutation(perm, electronsNumber), samples_);
                    if (settings.sortRemainder()){
                        auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                                subGroup->representative()->maximum().tail(electronsNumber - *countIterator).positionsVector(),
                                sortedGroup->representative()->maximum().tail(electronsNumber - *countIteratorSuperGroup).positionsVector());
                        subGroup->permuteAll(BestMatch::tailToFullPermutation(perm, electronsNumber), samples_);
                    }
                    sortedGroup->emplace_back(*subGroup);
                    isSimilarQ = true;
                    break;
                }
            }
            countIteratorSuperGroup++;
        }
        if (!isSimilarQ) {
            // adds subGroup as a new sortedGroup to superGroup
            superGroup.emplace_back(Group({*subGroup}));
            countsSuperGroup.emplace_back(*countIterator);
        }
        countIteratorSuperGroup = countsSuperGroup.begin();
        countIterator++;
    }
    group = superGroup;

    // sort by function value before leaving
    group.sort();
}

std::list<long> ReferencePositionsClusterer::getRelevantIndices(const ElectronsVector &electrons) {
    return NearestElectrons::getNearestElectronsIndices(electrons, nuclei_, positions_,
                                                        settings.maximalCount(), settings.maximalDistance(),
                                                        distanceFunction_, settings.valenceOnly());
}
