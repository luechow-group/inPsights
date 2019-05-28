//
// Created by heuer on 12.12.18.
//

#include <BestMatchDistanceIdentityClusterer.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <BestMatchDistance.h>
#include <ValueSorter.h>

namespace Settings {
    BestMatchDistanceSimilarityClusterer::BestMatchDistanceSimilarityClusterer()
    : ISettings(VARNAME(BestMatchDistanceSimilarityClusterer)) {
        similarityRadius.onChange_.connect(
                [&](double value) {
                    if(value < ::BestMatchDistanceIdentityClusterer::settings.identityRadius())
                        throw std::invalid_argument(
                                "The " + similarityRadius.name() + " with " + std::to_string(similarityRadius())
                                + " is smaller than the "+ ::BestMatchDistanceIdentityClusterer::settings.identityRadius.name() 
                                + " with "
                                + std::to_string(::BestMatchDistanceIdentityClusterer::settings.identityRadius()));
                });
    }

    BestMatchDistanceSimilarityClusterer::BestMatchDistanceSimilarityClusterer(const YAML::Node &node)
            : BestMatchDistanceSimilarityClusterer() {
        doubleProperty::decode(node, similarityRadius);
        doubleProperty::decode(node, similarityValueIncrement);
    }

    void BestMatchDistanceSimilarityClusterer::appendToNode(YAML::Node &node) const {
        node[className][similarityRadius.name()] = similarityRadius();
        node[className][similarityValueIncrement.name()] = similarityValueIncrement();
    }
}
YAML_SETTINGS_DEFINITION(Settings::BestMatchDistanceSimilarityClusterer)

Settings::BestMatchDistanceSimilarityClusterer BestMatchDistanceSimilarityClusterer::settings = Settings::BestMatchDistanceSimilarityClusterer();


BestMatchDistanceSimilarityClusterer::BestMatchDistanceSimilarityClusterer(std::vector<Sample> &samples)
        : samples_(samples){}
        
// assumes a sorted reference vector
void BestMatchDistanceSimilarityClusterer::cluster(Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    auto similarityRadius = settings.similarityRadius();
    auto valueIncrement = settings.similarityValueIncrement();


    // first, make sure group is sorted
    group.sort();

    // insert first element
    Group supergroup({Group({*group.begin()})});
    group.erase(group.begin());

    //Presort
    for (auto subgroup = group.begin(); subgroup != group.end(); ++subgroup) {

        // iterate over all supergroup members within the value range
        std::list<bool> outsideQ;
        for(auto subgroupFromSuperGroup = supergroup.begin();
        subgroupFromSuperGroup != supergroup.end(); ++subgroupFromSuperGroup) {

            auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
                    subgroup->representative()->maximum().positionsVector(),
                    subgroupFromSuperGroup->representative()->maximum().positionsVector());
            if(norm > similarityRadius)
                outsideQ.emplace_back(true);
            else
                outsideQ.emplace_back(false);
        }
        if(std::all_of(outsideQ.begin(), outsideQ.end(), [](bool b){return b;})) {
            supergroup.emplace_back(Group({*subgroup}));
            group.erase(subgroup);
            subgroup -= 1;
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
}
