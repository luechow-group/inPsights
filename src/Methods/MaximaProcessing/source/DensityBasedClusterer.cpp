//
// Created by heuer on 12.12.18.
// Edited by reuter on 21.05.18.
//

#include <DensityBasedClusterer.h>
#include <DistanceClusterer.h>
#include <BestMatchDistance.h>
#include <Enumerate.h>
#include <Reference.h>

namespace Settings {
    DensityBasedClusterer::DensityBasedClusterer()
    : ISettings(VARNAME(DensityBasedClusterer)) {
        radius.onChange_.connect(
                [&](double value) {
                    if(value < ::DistanceClusterer::settings.radius())
                        throw std::invalid_argument(
                                "The " + radius.name() + "=" + std::to_string(radius())
                                + " is smaller than the " + ::DistanceClusterer::settings.name()
                                + "::"  + ::DistanceClusterer::settings.radius.name()
                                + " of " + std::to_string(::DistanceClusterer::settings.radius()) + ".");
                });
        minimalClusterSize.onChange_.connect(
                [&](size_t value) {
                    if(!value >= 1)
                        throw std::invalid_argument(
                                "The " + minimalClusterSize.name() + "=" + std::to_string(minimalClusterSize())
                                + " must be 1 or greater.");
                });
    }

    DensityBasedClusterer::DensityBasedClusterer(const YAML::Node &node)
            : DensityBasedClusterer() {
        doubleProperty::decode(node, radius);
        size_tProperty ::decode(node, minimalClusterSize);
    }

    void DensityBasedClusterer::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius();
        node[className][minimalClusterSize.name()] = minimalClusterSize();
    }
}
YAML_SETTINGS_DEFINITION(Settings::DensityBasedClusterer)

Settings::DensityBasedClusterer DensityBasedClusterer::settings = Settings::DensityBasedClusterer();


DensityBasedClusterer::DensityBasedClusterer(std::vector<Sample> &samples)
        : samples_(samples) {};

double DensityBasedClusterer::wrapper(const Group &g1, const Group &g2) {
    return BestMatch::Distance::compare<Eigen::Infinity, 2>(
            g1.representative()->maximum().positionsVector(),
            g2.representative()->maximum().positionsVector()).metric;
};
        
void DensityBasedClusterer::cluster(Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    group.sortAll();

    auto eps = settings.radius();
    auto minPts = settings.minimalClusterSize();

    DensityBasedScan<double, Group, DensityBasedClusterer::wrapper> dbscan(group);
    auto result = dbscan.findClusters(eps, minPts);

    Group supergroup(static_cast<Group::size_type>(result.numberOfClusters));

    for (int i = 0; i < result.numberOfClusters; ++i)
        for (auto  [j, g] : enumerate(group))
            if (result.labels[j] == i)
                supergroup[i].emplace_back(std::move(g));

    orderByBestMatchDistance(supergroup, eps);

    group = supergroup;
}

void DensityBasedClusterer::orderByBestMatchDistance(Group &supergroup, double threshold) const {
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
                        j = subgroup.erase(j);

                        // the iterator has to be set back by one because the j element was erased and
                        // ++j of the for loop would otherwise skip one group of subgroup
                        --j;
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
