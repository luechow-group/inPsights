//
// Created by leonard on 13.05.19.
//

#include <LocalBondSimilarityClusterer.h>
#include <BestMatchDistance.h>
#include <BestMatch.h>
#include <NearestElectrons.h>
#include <Reference.h>
#include <Group.h>

namespace Settings {
    LocalBondSimilarityClusterer::LocalBondSimilarityClusterer()
            : ISettings(VARNAME(LocalBondSimilarityClusterer)) {};

    LocalBondSimilarityClusterer::LocalBondSimilarityClusterer(const YAML::Node &node)
            : LocalBondSimilarityClusterer() {
        doubleProperty::decode(node, similarityRadius);
        intProperty::decode(node, index1);
        intProperty::decode(node, index2);
        longProperty::decode(node, maximalCount);
    };

    void LocalBondSimilarityClusterer::appendToNode(YAML::Node &node) const {
        node[className][similarityRadius.name()] = similarityRadius();
        node[className][index1.name()] = index1();
        node[className][index2.name()] = index2();
        node[className][maximalCount.name()] = maximalCount();
    };
}

YAML_SETTINGS_DEFINITION(Settings::LocalBondSimilarityClusterer)

Settings::LocalBondSimilarityClusterer LocalBondSimilarityClusterer::settings = Settings::LocalBondSimilarityClusterer();

LocalBondSimilarityClusterer::LocalBondSimilarityClusterer(std::vector<Sample> &samples, AtomsVector nuclei)
        : samples_(samples),
          nuclei_(nuclei) {}

void LocalBondSimilarityClusterer::cluster(Group &group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto similarityRadius = settings.similarityRadius();
    long electronsNumber = group.representative()->maximum().numberOfEntities();

    group.sort();

    std::list<long> subIndices;
    Eigen::PermutationMatrix<Eigen::Dynamic> permutation;

    // sorting relevant electrons to the front
    for (auto subGroup = group.begin(); subGroup != group.end(); ++subGroup) {
        subIndices = LocalBondSimilarityClusterer::getRelevantIndices(subGroup->representative()->maximum());
        permutation = BestMatch::getPermutationToFront(subIndices, electronsNumber);
        subGroup->permuteAll(permutation, samples_);
    }

    Group superGroup({Group({*group.begin()})});

    bool isSimilarQ;

    for (auto subGroup = std::next(group.begin()); subGroup != group.end(); ++subGroup) {
        isSimilarQ = false;
        for (auto sortedGroup = superGroup.begin(); sortedGroup != superGroup.end(); ++sortedGroup) {
            auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subGroup->representative()->maximum().getFirstElements(settings.maximalCount()).positionsVector(),
                    sortedGroup->representative()->maximum().getFirstElements(settings.maximalCount()).positionsVector());

            if (norm < similarityRadius) {
                subGroup->permuteAll(BestMatch::getFullPermutation(perm, electronsNumber), samples_);
                sortedGroup->emplace_back(*subGroup);
                isSimilarQ = true;
                break;
            }
        }
        if (!isSimilarQ) {
            superGroup.emplace_back(Group({*subGroup}));
        }
    }
    group = superGroup;
}

std::list<long> LocalBondSimilarityClusterer::getRelevantIndices(const ElectronsVector &electrons) {
    Eigen::Vector3d position =
            (nuclei_[settings.index1()].position() + nuclei_[settings.index2()].position()) / 2;
    return NearestElectrons::getNearestValenceIndices(electrons, nuclei_, position, settings.maximalCount());
}
