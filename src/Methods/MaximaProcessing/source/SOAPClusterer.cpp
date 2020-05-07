/* Copyright (C) 2018-2019 Michael Heuer.
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

#include "SOAPClusterer.h"
#include "PreClusterer.h"
#include <spdlog/spdlog.h>
#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>
#include <Reference.h>

using namespace SOAP;

namespace Settings {
    SOAPClusterer::SOAPClusterer()
            : ISettings(VARNAME(SOAPClusterer)) {
        similarityThreshold.onChange_.connect(
                [&](double value) {
                    if(value <= 0 || value > 1)
                        throw std::invalid_argument(
                                "The " + similarityThreshold.name() + " with " + std::to_string(similarityThreshold())
                                + " must be within the range (0,1].");
                });
    }

    SOAPClusterer::SOAPClusterer(const YAML::Node &node)
            : SOAPClusterer() {
        doubleProperty::decode(node, similarityThreshold);
        doubleProperty::decode(node, distanceMatrixCovarianceTolerance);
    }

    void SOAPClusterer::appendToNode(YAML::Node &node) const {
        node[className][similarityThreshold.name()] = similarityThreshold();
        node[className][distanceMatrixCovarianceTolerance.name()] = distanceMatrixCovarianceTolerance();
    }
}
YAML_SETTINGS_DEFINITION(Settings::SOAPClusterer)

Settings::SOAPClusterer SOAPClusterer::settings = Settings::SOAPClusterer();


SOAPClusterer::SOAPClusterer(const AtomsVector& atoms, std::vector<Sample> &samples)
    : IClusterer(samples), atoms_(atoms) {
    ParticleKit::create(atoms_, (*samples.begin()).sample_);
};

void SOAPClusterer::cluster(Group& group){
    assert(!group.isLeaf() && "The group cannot be a leaf.");

    // TODO use function value to narrow down guesses => presort lists

    // Calculate spectra
    spdlog::info("Calculating {} spectra in {} mode...",
            group.size(),SOAP::General::toString(SOAP::General::settings.mode.get()));
#pragma omp parallel for default(none) shared(atoms_, group)
    for (auto it = group.begin(); it < group.end(); ++it) {

        it->representative()->setSpectrum(MolecularSpectrum(
                {atoms_, it->representative()->maximum()}));
        spdlog::info("calculated spectrum {}", std::distance(group.begin(), it));
    }

    auto similarityThreshold = settings.similarityThreshold();
    auto toleranceRadius = settings.distanceMatrixCovarianceTolerance();
    auto numericalPrecisionEpsilon = SOAP::General::settings.comparisonEpsilon.get();

    spdlog::debug("Group before start: {}", ToString::groupToString(group));

    Group supergroup({{*group.begin()}});
    spdlog::debug("Supergroup before start: {}", ToString::groupToString(supergroup));

    for(auto groupIt = std::next(group.begin()); groupIt !=group.end(); groupIt++) {
        spdlog::info("{} out of {}",
                std::distance(group.begin(), groupIt), std::distance(group.begin(), group.end()));

        bool foundMatchQ = false;

        spdlog::debug("  Outer loop groupIt {}: {}",
                std::distance(group.begin(), groupIt), ToString::groupToString(*groupIt));
        // check if current group matches any of the supergroup subgroups
        for(auto subgroupOfSupergroupIt = supergroup.begin(); subgroupOfSupergroupIt != supergroup.end(); ++subgroupOfSupergroupIt) {

            spdlog::debug("    Inner loop subgroupOfSupergroupIt {}: {}",
                    std::distance(supergroup.begin(),subgroupOfSupergroupIt), ToString::groupToString(*subgroupOfSupergroupIt));

            assert(!groupIt->representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");
            assert(!subgroupOfSupergroupIt->representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");

            auto comparisionResult = BestMatch::SOAPSimilarity::compare(
                    groupIt->representative()->spectrum(),
                    subgroupOfSupergroupIt->representative()->spectrum(),
                    toleranceRadius,
                    similarityThreshold, numericalPrecisionEpsilon);

            spdlog::info("  comparing it with {} out of {}: {}",
                    std::distance(supergroup.begin(), subgroupOfSupergroupIt)+1,
                    std::distance(supergroup.begin(), supergroup.end()),
                    comparisionResult.metric);

            // if so, put permute the current group and put it into the supergroup subgroup and stop searching
            if (comparisionResult.metric >= (similarityThreshold-numericalPrecisionEpsilon)) {
                groupIt->permuteAll(comparisionResult.permutation, samples_);

                subgroupOfSupergroupIt->emplace_back(*groupIt);

                spdlog::debug("    Match: Inner loop subgroupOfSupergroupIt {}: {}",
                              std::distance(supergroup.begin(),subgroupOfSupergroupIt),
                              ToString::groupToString(*subgroupOfSupergroupIt));

                foundMatchQ = true;
                break;
            }
            spdlog::debug("    No match.");

        }
        if(!foundMatchQ) {
            supergroup.emplace_back(Group({*groupIt}));

            spdlog::debug(" No match found. Group {} was {} to supergroup: {}",
                          std::distance(group.begin(), groupIt),
                          ToString::groupToString(*groupIt),
                          ToString::groupToString(supergroup));
        }

    }
    group = supergroup;

    spdlog::debug("Result after loop: {}", ToString::groupToString(group));

    // sort by function value before leaving
    group.sortAll();

    spdlog::debug("Final result after sorting: {}", ToString::groupToString(group));
}