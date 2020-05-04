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
#include <StructuralSimilarity.h>
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

    Group supergroup({{*group.begin()}});

    for(auto groupIt = std::next(group.begin()); groupIt !=group.end(); groupIt++) {
        spdlog::info("{} out of {}", std::distance(group.begin(), groupIt), std::distance(group.begin(), group.end()));
        bool foundMatchQ = false;

        // check if current group matches any of the supergroup subgroups
        for(auto subgroupOfSupergroupIt = supergroup.begin(); subgroupOfSupergroupIt != supergroup.end(); ++subgroupOfSupergroupIt) {

            assert(!groupIt->representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");
            assert(!subgroupOfSupergroupIt->representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");

            auto comparisionResult = BestMatch::SOAPSimilarity::compare(
                    groupIt->representative()->spectrum(),
                    subgroupOfSupergroupIt->representative()->spectrum(),
                    toleranceRadius,
                    similarityThreshold);

            spdlog::info("  comparing it with {} out of {}: {}",
                    std::distance(supergroup.begin(), subgroupOfSupergroupIt),
                    std::distance(supergroup.begin(), supergroup.end()),
                    comparisionResult.metric);

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

    // sort by function value before leaving
    group.sortAll();
}