// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
        maximalStructuresNumber.onChange_.connect(
                [&](unsigned value) {
                    if (value < 1)
                        throw std::invalid_argument("The maximal number of structures to print cannot be less than 1.");
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
        unsignedProperty::decode(node[className], maximalStructuresNumber);
        doubleProperty ::decode(node[className], motifThreshold);
        boolProperty ::decode(node[className], deleteCoreElectrons);
        boolProperty ::decode(node[className], printAllMaxima);
        boolProperty ::decode(node[className], doEnergyPartitioning);
        boolProperty ::decode(node[className], calculateSpinCorrelations);
        boolProperty ::decode(node[className], shuffleMaxima);
        if (doEnergyPartitioning.get() and not calculateSpinCorrelations.get())
            throw std::invalid_argument("calcSpinCorrelation has to be true, if doEnergyPartitioning is true.");
    }

    void MaximaProcessing::appendToNode(YAML::Node &node) const {
        node[className][seed.name()] = seed.get();
        node[className][binaryFileBasename.name()] = binaryFileBasename.get();
        node[className][samplesToAnalyze.name()] = samplesToAnalyze.get();
        node[className][minimalClusterWeight.name()] = minimalClusterWeight.get();
        node[className][maximalStructuresNumber.name()] = maximalStructuresNumber.get();
        node[className][motifThreshold.name()] = motifThreshold.get();
        node[className][deleteCoreElectrons.name()] = deleteCoreElectrons.get();
        node[className][printAllMaxima.name()] = printAllMaxima.get();
        node[className][doEnergyPartitioning.name()] = doEnergyPartitioning.get();
        node[className][calculateSpinCorrelations.name()] = calculateSpinCorrelations.get();
        node[className][shuffleMaxima.name()] = shuffleMaxima.get();
    }
}

YAML_SETTINGS_DEFINITION(Settings::MaximaProcessing)

namespace MaximaProcessing {
    Settings::MaximaProcessing settings = Settings::MaximaProcessing();
}