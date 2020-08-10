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

#include "MaximaProcessingSettings.h"
#include <random>
#include <spdlog/spdlog.h>

namespace Settings {
    MaximaProcessing::MaximaProcessing()
    : ISettings(VARNAME(MaximaProcessing)) {
        seed.onChange_.connect(
                [&](unsigned value) {
                    if (value == 0) {
                        seed = std::random_device()();
                        spdlog::info("Generated seed from random device.");
                    }
                    spdlog::info("Seed: {}", seed());
                });
        samplesToAnalyze.onChange_.connect(
                [&](unsigned value) {
                    if (not (value >= 0 and value < std::numeric_limits<unsigned>::max()))
                        throw std::invalid_argument("The samples to analyze must not be negative.");
                });
        minimalClusterWeight.onChange_.connect(
                [&](double value) {
                    if (value < 0.0)
                        throw std::invalid_argument("The minimal cluster weight cannot be negative.");
                    else if (value >= 1.0)
                        throw std::invalid_argument("The minimal cluster weight cannot be 1 or larger.");
                });
        maximalNumberOfStructuresToPrint.onChange_.connect(
                [&](unsigned value) {
                    if (value < 0)
                        throw std::invalid_argument("The maximal number of structures to print cannot be negative.");
                });

        motifThreshold.onChange_.connect(
                [&](double value) {
                    if (value < 0.0)
                        throw std::invalid_argument("The motif threshold cannot be negative.");
                    else if (value > 1.0)
                        throw std::invalid_argument("The motif threshold cannot greater than 1.");
                });
    }

    MaximaProcessing::MaximaProcessing(const YAML::Node &node)
    : MaximaProcessing() {
        YAML::convert<Property<std::string>>::decode(node[className], binaryFileBasename);
        unsignedProperty::decode(node[className], seed);
        unsignedProperty::decode(node[className], samplesToAnalyze);
        doubleProperty ::decode(node[className], minimalClusterWeight);
        unsignedProperty::decode(node[className], maximalNumberOfStructuresToPrint);
        doubleProperty ::decode(node[className], motifThreshold);
        boolProperty ::decode(node[className], deleteCoreElectrons);
    }

    void MaximaProcessing::appendToNode(YAML::Node &node) const {
        node[className][seed.name()] = seed.get();
        node[className][binaryFileBasename.name()] = binaryFileBasename.get();
        node[className][samplesToAnalyze.name()] = samplesToAnalyze.get();
        node[className][minimalClusterWeight.name()] = minimalClusterWeight.get();
        node[className][maximalNumberOfStructuresToPrint.name()] = maximalNumberOfStructuresToPrint.get();
        node[className][motifThreshold.name()] = motifThreshold.get();
        node[className][deleteCoreElectrons.name()] = deleteCoreElectrons.get();
    }
}

YAML_SETTINGS_DEFINITION(Settings::MaximaProcessing)

namespace MaximaProcessing {
    Settings::MaximaProcessing settings = Settings::MaximaProcessing();
}