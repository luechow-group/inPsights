//
// Created by heuer on 12.12.18.
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
                    if(value > ::BestMatchDistanceSimilarityClusterer::settings.similarityRadius())
                        throw std::invalid_argument(
                                "The " + clusterRadius.name() + " with " + std::to_string(clusterRadius())
                                + " is greater than the "+ ::BestMatchDistanceSimilarityClusterer::settings.similarityRadius.name() 
                                + " with "
                                + std::to_string(::BestMatchDistanceSimilarityClusterer::settings.similarityRadius()));
                });
    }

    BestMatchDistanceDensityBasedClusterer::BestMatchDistanceDensityBasedClusterer(const YAML::Node &node)
            : BestMatchDistanceDensityBasedClusterer() {
        doubleProperty::decode(node[className], clusterRadius);
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
            g1.representative()->maximum(),
            g2.representative()->maximum()).metric;
};

bool BestMatchDistanceDensityBasedClusterer::SortElement::operator<(const SortElement &rhs) const {
    return bestMatch_.metric < rhs.bestMatch_.metric;
}


BestMatchDistanceDensityBasedClusterer::SortElement::SortElement(BestMatch::Result bestMatch, Group::iterator it)
        : bestMatch_(std::move(bestMatch)), it_(it) {}
        
void BestMatchDistanceDensityBasedClusterer::cluster(Group& group) {
    group.sortAll();
    auto threshold = settings.clusterRadius() * 2 + 0.01; // TODO WHY IS THIS CORRECTION NECESSARY?

    DensityBasedScan<double, Group, BestMatchDistanceDensityBasedClusterer::wrapper> dbscan(group);
    auto result = dbscan.findClusters(threshold, 1); // why multiplication by 2 is needed?

    Group supergroup(static_cast<Group::size_type>(result.numberOfClusters));

    for (int i = 0; i < result.numberOfClusters; ++i)
        for (auto  [j, g] : enumerate(group))
            if (result.labels[j] == i)
                supergroup[i].emplace_back(std::move(g));

    orderByBestMatchDistance(supergroup);

    group = supergroup;
}

// TODO this is a more general group sorting method
void BestMatchDistanceDensityBasedClusterer::orderByBestMatchDistance(Group &supergroup) const {// order clusters by best match distance //TODO make specialized sorting method
    for (auto &subgroup : supergroup) {
        sort(subgroup.begin(), subgroup.end());

        // iterate over all similarReferences in the cluster
        for (auto i = subgroup.begin(); i != subgroup.end(); ++i) {
            std::vector<SortElement> bestMatchDistances;

            // make a list of best match distances and permutations starting with the next similarReferences object
            for (auto j = i + 1; j != subgroup.end(); ++j) {
                bestMatchDistances.emplace_back(
                        SortElement(BestMatch::Distance::compare<Eigen::Infinity, 2>(
                                j.base()->representative()->maximum(),
                                i.base()->representative()->maximum()), j)
                );
            }

            // check if the list contains more than one element
            if (bestMatchDistances.size() > 1) {
                // find the SimilarReferences object whose representativeReference is closest to the i
                auto minIt = min_element(bestMatchDistances.begin(), bestMatchDistances.end());

                // permute and swap
                minIt.base()->it_.base()->permuteAll(minIt.base()->bestMatch_.permutation, samples_);
                if (i + 1 != minIt.base()->it_) {
                    iter_swap(i + 1, minIt.base()->it_);
                }
            } else if (bestMatchDistances.size() == 1) { // only one element left
                (i + 1).base()->permuteAll(bestMatchDistances[0].bestMatch_.permutation, samples_);
            }
        }
    }
}
