//
// Created by leonard on 13.05.19.
//

#include <LocalBondSimilarityClusterer.h>
#include <BestMatchDistance.h>
#include <BestMatch.h>
#include <NearestElectrons.h>
#include <Reference.h>
#include <Group.h>
#include <functional>

namespace Settings {
    LocalBondSimilarityClusterer::LocalBondSimilarityClusterer()
            : ISettings(VARNAME(LocalBondSimilarityClusterer)) {
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
        maximalCount.onChange_.connect(
                [&](long value) {
                    if (not (value > 0))
                        throw std::invalid_argument("The maximalCount has to be larger than zero.");
                });
    };

    LocalBondSimilarityClusterer::LocalBondSimilarityClusterer(const YAML::Node &node)
            : LocalBondSimilarityClusterer() {
        doubleProperty::decode(node, similarityRadius);
        longProperty::decode(node, maximalCount);
        stringProperty::decode(node, distanceMode);
    };

    void LocalBondSimilarityClusterer::appendToNode(YAML::Node &node) const {
        node[className][similarityRadius.name()] = similarityRadius();
        node[className][maximalCount.name()] = maximalCount();
        node[className][distanceMode.name()] = distanceMode();
    };
}

YAML_SETTINGS_DEFINITION(Settings::LocalBondSimilarityClusterer)

Settings::LocalBondSimilarityClusterer LocalBondSimilarityClusterer::settings = Settings::LocalBondSimilarityClusterer();

LocalBondSimilarityClusterer::LocalBondSimilarityClusterer(std::vector<Sample> &samples, AtomsVector &nuclei,
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
    return NearestElectrons::getNearestValenceIndices(electrons, nuclei_, positions_,
            settings.maximalCount(), distanceFunction_);
}
