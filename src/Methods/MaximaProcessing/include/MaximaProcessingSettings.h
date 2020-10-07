// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MAXIMAPROCESSINGSETTINGS_H
#define INPSIGHTS_MAXIMAPROCESSINGSETTINGS_H

#include <ISettings.h>
#include <Property.h>
#include <string>

enum class SamplesToAnalyze{
    small = 10000,
    medium = 100000,
    large = 1000000,
    all = 0 };

namespace Settings {
    class MaximaProcessing : public ISettings {
    public:
        Property<unsigned> seed = {12345, VARNAME(seed)}; // A random seed is picked for a value of 0
        Property<std::string> binaryFileBasename = {"(name of the .yml input file)", VARNAME(binaryFileBasename)};
        Property<unsigned> samplesToAnalyze = {unsigned(SamplesToAnalyze::small), VARNAME(samplesToAnalyze)};
        Property<double> minimalClusterWeight = {0.01, VARNAME(minimalClusterWeight)};
        Property<unsigned> maximalNumberOfStructuresToPrint = {std::numeric_limits<unsigned>::max(), VARNAME(maximalNumberOfStructuresToPrint)};
        Property<double> motifThreshold = {1.00, VARNAME(motifThreshold)};
        Property<bool> deleteCoreElectrons = {false, VARNAME(deleteCoreElectrons)};
        Property<bool> printAllMaxima = {false, VARNAME(printAllMaxima)};

        MaximaProcessing();
        explicit MaximaProcessing(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::MaximaProcessing)

namespace MaximaProcessing {
    extern Settings::MaximaProcessing settings;
}

#endif //INPSIGHTS_MAXIMAPROCESSINGSETTINGS_H
