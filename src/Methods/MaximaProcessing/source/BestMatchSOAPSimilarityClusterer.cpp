//
// Created by heuer on 12.12.18.
//

#include "BestMatchSOAPSimilarityClusterer.h"
#include "BestMatchDistanceSimilarityClusterer.h"
#include <StructuralSimilarity.h>
#include <spdlog/spdlog.h>
#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>
#include <Reference.h>

using namespace SOAP;

namespace Settings {
    BestMatchSOAPSimilarityClusterer::BestMatchSOAPSimilarityClusterer()
            : ISettings(VARNAME(BestMatchSOAPSimilarityClusterer)) {
        threshold.onChange_.connect(
                [&](double value) {
                    if(value <= 0 || value > 1)
                        throw std::invalid_argument(
                                "The " + threshold.name() + " with " + std::to_string(threshold())
                                + " must be within the range (0,1].");
                });
    }

    BestMatchSOAPSimilarityClusterer::BestMatchSOAPSimilarityClusterer(const YAML::Node &node)
            : BestMatchSOAPSimilarityClusterer() {
        doubleProperty::decode(node, threshold);
    }

    void BestMatchSOAPSimilarityClusterer::appendToNode(YAML::Node &node) const {
        node[className][threshold.name()] = threshold();
    }
}
YAML_SETTINGS_DEFINITION(Settings::BestMatchSOAPSimilarityClusterer)

Settings::BestMatchSOAPSimilarityClusterer BestMatchSOAPSimilarityClusterer::settings = Settings::BestMatchSOAPSimilarityClusterer();


BestMatchSOAPSimilarityClusterer::BestMatchSOAPSimilarityClusterer(
        const AtomsVector& atoms,
        std::vector<Sample> &samples)
        : atoms_(atoms), samples_(samples) {
    ParticleKit::create(atoms_, (*samples.begin()).sample_);
};

void BestMatchSOAPSimilarityClusterer::cluster(Group& group){
    assert(!group.isLeaf() && "The group cannot be a leaf.");

    group.sortAll(); // sort, to make sure that the most probable structures are the representatives

    // Calculate spectra
    spdlog::info("Calculating {} spectra...", group.size());
#pragma omp parallel for default(none) shared(atoms_, group)
    for (auto it = group.begin(); it < group.end(); ++it) {
        it->representative()->setSpectrum(MolecularSpectrum({atoms_, it->representative()->maximum()}));
        spdlog::info("calculated spectrum {}", std::distance(group.begin(), it));
    }

    auto similarityThreshold = settings.threshold();

    Group supergroup({{*group.begin()}});

    for(auto groupIt = std::next(group.begin()); groupIt !=group.end(); groupIt++) {
        spdlog::info("{} out of {}", std::distance(group.begin(), groupIt), std::distance(group.begin(), group.end()));
        bool foundMatchQ = false;

        // check if current group matches any of the supergroup subgroups
        for(auto subgroupOfSupergroupIt = supergroup.begin(); subgroupOfSupergroupIt != supergroup.end(); ++subgroupOfSupergroupIt) {
            spdlog::info("  comparing it with {} out of {}",
                    std::distance(supergroup.begin(), subgroupOfSupergroupIt),
                    std::distance(supergroup.begin(), supergroup.end()));
            assert(!groupIt->representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");
            assert(!subgroupOfSupergroupIt->representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");

            auto comparisionResult = BestMatch::SOAPSimilarity::compare(
                    groupIt->representative()->spectrum(),
                    subgroupOfSupergroupIt->representative()->spectrum(),
                    BestMatchDistanceSimilarityClusterer::settings.similarityRadius(),//TODO  CAREFUL!!! ADD OWN PARAMETER FOR TOLSIM RADIUS
                    settings.threshold());

            // if so, put permute the current group and put it into the supergroup subgroup and stop searching
            if (comparisionResult.metric >= similarityThreshold) {
                groupIt->permuteAll(comparisionResult.permutation, samples_);

                subgroupOfSupergroupIt->emplace_back(*groupIt);

                foundMatchQ = true;
                break;
            }
        }
        if(!foundMatchQ)
            supergroup.emplace_back(Group({*groupIt}));
    }
    group = supergroup;
}