// Copyright (C) 2018-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
        doubleProperty::decode(node, maxValueDelta);

    }

    void SOAPClusterer::appendToNode(YAML::Node &node) const {
        node[className][similarityThreshold.name()] = similarityThreshold();
        node[className][distanceMatrixCovarianceTolerance.name()] = distanceMatrixCovarianceTolerance();
        node[className][maxValueDelta.name()] = maxValueDelta();
    }
}
YAML_SETTINGS_DEFINITION(Settings::SOAPClusterer)

Settings::SOAPClusterer SOAPClusterer::settings = Settings::SOAPClusterer();


SOAPClusterer::SOAPClusterer(const AtomsVector& atoms, std::vector<Sample> &samples)
    : IClusterer(samples), atoms_(atoms) {
    ParticleKit::create(atoms_, (*samples.begin()).sample_);
};

void SOAPClusterer::cluster(Cluster& cluster){
    assert(!cluster.isLeaf() && "The cluster cannot be a leaf.");

    // TODO use function value to narrow down guesses => presort lists

    // Calculate spectra
    spdlog::info("Calculating {} spectra in {} mode...",
            cluster.size(),SOAP::General::toString(SOAP::General::settings.mode.get()));
#pragma omp parallel for default(none) shared(atoms_, cluster)
    for (auto it = cluster.begin(); it < cluster.end(); ++it) {
        it->representative()->setSpectrum(MolecularSpectrum(
                {atoms_, it->representative()->maximum()}));
        spdlog::info("calculated spectrum {}", std::distance(cluster.begin(), it));
    }

    auto similarityThreshold = settings.similarityThreshold();
    auto toleranceRadius = settings.distanceMatrixCovarianceTolerance();
    auto numericalPrecisionEpsilon = SOAP::General::settings.comparisonEpsilon();
    auto maxValueDelta = SOAPClusterer::settings.maxValueDelta();

    spdlog::debug("Cluster before start: {}", ToString::clusterToString(cluster));

    Cluster supercluster;
    spdlog::debug("Supercluster before start: {}", ToString::clusterToString(supercluster));
    for(const auto& [i, subcluster] : enumerate(cluster)) {
        spdlog::info("{} out of {}", i+1, cluster.size());

        bool foundMatchQ = false;

        spdlog::debug("  Outer loop clusterIt {}: {}", i, ToString::clusterToString(subcluster));
        // check if current cluster matches any of the supercluster subclusters
        for(const auto & [j, subclusterOfSupercluster] : enumerate(supercluster)){
            spdlog::debug("    Inner loop subclusterOfSuperclusterIt {}: {}", j, ToString::clusterToString(subclusterOfSupercluster));

            assert(!subcluster.representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");
            assert(!subclusterOfSupercluster.representative()->spectrum().molecularCenters_.empty() && "Spectrum cannot be empty.");


            spdlog::debug("    Supercluster status before comparison: {}", ToString::clusterToString(supercluster));

            if( std::abs(subclusterOfSupercluster.representative()->value() - subcluster.representative()->value()) < maxValueDelta) {
                auto comparisionResult = BestMatch::SOAPSimilarity::compare(
                        subcluster.representative()->spectrum(),
                        subclusterOfSupercluster.representative()->spectrum(),
                        toleranceRadius, similarityThreshold, numericalPrecisionEpsilon);

                spdlog::debug("    Supercluster status after comparison: {}", ToString::clusterToString(supercluster));

                spdlog::info("  comparing it with {} out of {}: {}",
                             j + 1, supercluster.size(),
                             comparisionResult.metric);

                // if so, put permute the current cluster and put it into the supercluster subcluster and stop searching
                if (comparisionResult.metric >= (similarityThreshold - numericalPrecisionEpsilon)) {
                    auto permutations = subcluster.representative()->spectrum().molecule_.splitAllParticlePermutation(comparisionResult.permutation);
                    subcluster.permuteAll(permutations, samples_);

                    supercluster[j].emplace_back(subcluster);

                    spdlog::debug("    Match: Inner loop subclusterOfSuperclusterIt {}: {}",
                                  j + 1, ToString::clusterToString(subclusterOfSupercluster));
                    spdlog::debug("    Match. End of inner loop. Supercluster status: {}",
                                  ToString::clusterToString(supercluster));
                    foundMatchQ = true;
                    break;
                }
            }

            spdlog::debug("    No match. End of inner loop. Supercluster status: {}", ToString::clusterToString(supercluster));
        }
        if(!foundMatchQ) {
            supercluster.emplace_back(Cluster({subcluster}));

            spdlog::debug(" No match found. Cluster {} {} was added to supercluster: {}",
                          i, ToString::clusterToString(subcluster), ToString::clusterToString(supercluster));
        }

    }
    cluster = supercluster;

    spdlog::debug("Result after loop: {}", ToString::clusterToString(cluster));

    // sort by function value before leaving
    cluster.sortAll();

    spdlog::debug("Final result after sorting: {}", ToString::clusterToString(cluster));
}
