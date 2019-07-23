//
// Created by heuer on 12.12.18.
//

#include <BestMatchDistanceIdentityClusterer.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <BestMatchDistance.h>
#include <ValueSorter.h>
#include <spdlog/spdlog.h>

namespace Settings {
    BestMatchDistanceIdentityClusterer::BestMatchDistanceIdentityClusterer()
    : ISettings(VARNAME(BestMatchDistanceIdentityClusterer)) {}

    BestMatchDistanceIdentityClusterer::BestMatchDistanceIdentityClusterer(const YAML::Node &node)
            : BestMatchDistanceIdentityClusterer() {
        doubleProperty::decode(node, identityRadius);
        doubleProperty::decode(node, identityValueIncrement);
    }

    void BestMatchDistanceIdentityClusterer::appendToNode(YAML::Node &node) const {
        node[className][identityRadius.name()] = identityRadius();
        node[className][identityValueIncrement.name()] = identityValueIncrement();
    }
}
YAML_SETTINGS_DEFINITION(Settings::BestMatchDistanceIdentityClusterer)

Settings::BestMatchDistanceIdentityClusterer BestMatchDistanceIdentityClusterer::settings = Settings::BestMatchDistanceIdentityClusterer();


BestMatchDistanceIdentityClusterer::BestMatchDistanceIdentityClusterer(std::vector<Sample> &samples)
        : samples_(samples) {}

void BestMatchDistanceIdentityClusterer::cluster(Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto identityRadius = settings.identityRadius();
    auto valueIncrement = settings.identityValueIncrement();

    // first, sort references by value
    group.sortAll();

    auto beginIt = group.begin();

    while (beginIt != group.end()) {
        auto total = group.size();//std::distance(group.begin(), group.end());
        auto endIt = std::upper_bound(beginIt, group.end(), Group(Reference(beginIt->representative()->value() + valueIncrement)));

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
}

void BestMatchDistanceIdentityClusterer::subLoop(Group& group,
        Group::iterator &beginIt,
        Group::iterator &it,
        Group::iterator &endIt,
        double distThresh,
        double valueIncrement) {

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
                    Group(Reference(beginIt->representative()->value() + valueIncrement)));
        } else it++;
    } else {  // don't consider spin flip
        if (norm <= distThresh) {
            addReference(group, beginIt, it, perm);
            endIt = std::upper_bound(beginIt, group.end(),
                    Group(Reference(beginIt->representative()->value() + valueIncrement)));
        } else it++;
    }
}

// TODO This method should be located inside of a reference container class
void BestMatchDistanceIdentityClusterer::addReference(Group& group,
        const Group::iterator &beginIt,
        Group::iterator &it,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const {

    it->permuteAll(bestMatch, samples_);
    beginIt->representative()->mergeReference(it);
    it = group.erase(it); // erase returns the iterator of the following element
}