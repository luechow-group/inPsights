//
// Created by heuer on 12.12.18.
// Edited by reuter on 21.05.18.
//

#include <BestMatchDistanceDensityBasedClusterer.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <BestMatchDistance.h>
#include <Enumerate.h>
#include <Reference.h>

namespace Settings {
    BestMatchDistanceDensityBasedClusterer::BestMatchDistanceDensityBasedClusterer()
    : ISettings(VARNAME(BestMatchDistanceDensityBasedClusterer)) {
        clusterRadius.onChange_.connect(
                [&](double value) {
                    if(value < ::BestMatchDistanceSimilarityClusterer::settings.similarityRadius())
                        throw std::invalid_argument(
                                "The " + clusterRadius.name() + " with " + std::to_string(clusterRadius())
                                + " is smaller than the "+ ::BestMatchDistanceSimilarityClusterer::settings.similarityRadius.name()
                                + " with "
                                + std::to_string(::BestMatchDistanceSimilarityClusterer::settings.similarityRadius()));
                });
    }

    BestMatchDistanceDensityBasedClusterer::BestMatchDistanceDensityBasedClusterer(const YAML::Node &node)
            : BestMatchDistanceDensityBasedClusterer() {
        doubleProperty::decode(node, clusterRadius);
    }

    void BestMatchDistanceDensityBasedClusterer::appendToNode(YAML::Node &node) const {
        node[className][clusterRadius.name()] = clusterRadius();
    }
}
YAML_SETTINGS_DEFINITION(Settings::BestMatchDistanceDensityBasedClusterer)

Settings::BestMatchDistanceDensityBasedClusterer BestMatchDistanceDensityBasedClusterer::settings = Settings::BestMatchDistanceDensityBasedClusterer();


BestMatchDistanceDensityBasedClusterer::BestMatchDistanceDensityBasedClusterer(std::vector<Sample> &samples)
        : samples_(samples) {};

double BestMatchDistanceDensityBasedClusterer::wrapper(const Group &g1, const Group &g2) {
    return BestMatch::Distance::compare<Eigen::Infinity, 2>(
            g1.representative()->maximum().positionsVector(),
            g2.representative()->maximum().positionsVector()).metric;
};
        
void BestMatchDistanceDensityBasedClusterer::cluster(Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    group.sortAll();

    auto threshold = settings.clusterRadius()*2;
    DensityBasedScan<double, Group, BestMatchDistanceDensityBasedClusterer::wrapper> dbscan(group);
    auto result = dbscan.findClusters(threshold, 1);// why multiplication by 2 is needed?

    Group supergroup(static_cast<Group::size_type>(result.numberOfClusters));

    for (int i = 0; i < result.numberOfClusters; ++i)
        for (auto  [j, g] : enumerate(group))
            if (result.labels[j] == i)
                supergroup[i].emplace_back(std::move(g));

    orderByBestMatchDistance(supergroup, threshold);

    group = supergroup;
}

void BestMatchDistanceDensityBasedClusterer::orderByBestMatchDistance(Group &supergroup, double threshold) const {
    for (auto &subgroup : supergroup) {
        sort(subgroup.begin(), subgroup.end());

        // starting sortedGroup with one active group, which is erased from subgroup
        Group sortedGroup({*subgroup.begin()});
        subgroup.erase(subgroup.begin());
        long activeGroups = 1;

        while (!subgroup.empty()){
            // setting newGroups empty again
            Group newGroups;

            // iterating over all active groups (at the end of sortedGroup)
            for (auto i = sortedGroup.end() - activeGroups; i != sortedGroup.end(); ++i){

                // iterating over all groups remaining in the unsorted subgroup
                for (auto j = subgroup.begin(); j != subgroup.end(); ++j) {
                    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                                    j->representative()->maximum().positionsVector(),
                                    i->representative()->maximum().positionsVector());

                    if (norm <= threshold) {
                        j->permuteAll(perm, samples_);

                        // moving j from subgroup to newGroups
                        newGroups.emplace_back(*j);
                        subgroup.erase(j);
                        j -= 1;
                    };
                };
            };

            // moving all groups from newGroups to sortedGroup()
            activeGroups = newGroups.size();
            for (auto &newGroup : newGroups) {
                sortedGroup.emplace_back(newGroup);
            };
        };

        subgroup = sortedGroup;
    }
}
