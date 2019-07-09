//
// Created by leonard on 13.05.19.
//

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
        similarityRadius.onChange_.connect(
                [&](double value) {
                    if (not (value > 0.0))
                        throw std::invalid_argument("The similarityRadius has to be larger than zero.");
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
        doubleProperty::decode(node, similarityRadius);
        doubleProperty::decode(node, maximalDistance);
        longProperty::decode(node, maximalCount);
        stringProperty::decode(node, distanceMode);
        boolProperty::decode(node, invertSelection);
        boolProperty::decode(node, valenceOnly);
    };

    void ReferencePositionsClusterer::appendToNode(YAML::Node &node) const {
        node[className][similarityRadius.name()] = similarityRadius();
        node[className][maximalDistance.name()] = maximalDistance();
        node[className][maximalCount.name()] = maximalCount();
        node[className][distanceMode.name()] = distanceMode();
        node[className][invertSelection.name()] = invertSelection();
        node[className][valenceOnly.name()] = valenceOnly();
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
        distanceFunction_ = Metrics::averageDistance;
    }
    else if (settings.distanceMode() == "minimum"){
        distanceFunction_ = Metrics::minimalDistance;
    }
}

void ReferencePositionsClusterer::cluster(Group &group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto similarityRadius = settings.similarityRadius();
    long electronsNumber = group.representative()->maximum().numberOfEntities();

    // sorts 'group' by the value ( -ln(|\Psi|^2) )
    group.sort();

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
                        subGroup->representative()->maximum().getFirstElements(*countIterator).positionsVector(),
                        sortedGroup->representative()->maximum().getFirstElements(*countIteratorSuperGroup).positionsVector());

                if (norm < similarityRadius) {
                    subGroup->permuteAll(BestMatch::getFullPermutation(perm, electronsNumber), samples_);
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
}

std::list<long> ReferencePositionsClusterer::getRelevantIndices(const ElectronsVector &electrons) {
    return NearestElectrons::getNearestElectronsIndices(electrons, nuclei_, positions_,
                                                        settings.maximalCount(), settings.maximalDistance(),
                                                        distanceFunction_, settings.valenceOnly());
}
