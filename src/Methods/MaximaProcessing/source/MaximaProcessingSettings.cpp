//
// Created by Michael Heuer on 2018-12-28.
//

#include "MaximaProcessingSettings.h"
#include <assert.h>
#include <spdlog/spdlog.h>

namespace Settings {
    MaximaProcessing::MaximaProcessing()
    : ISettings(VARNAME(MaximaProcessing)) {
        samplesToAnalyze.onChange_.connect(
                [&](unsigned value) {
                    if(value < std::numeric_limits<unsigned>::max())
                        spdlog::info("Analyzing {} samples.", value);
                    else if (value == 0)
                        spdlog::info("Analyzing all samples.");
                });
        minimalClusterWeight.onChange_.connect(
                [&](double value) {
                    if (value < 0.0)
                        throw std::invalid_argument("The minimal cluster weight cannot be negative.");
                    else if (value >= 1.0)
                        throw std::invalid_argument("The minimal cluster weight cannot be 1 or larger.");
                });
    }

    MaximaProcessing::MaximaProcessing(const YAML::Node &node)
    : MaximaProcessing() {
        YAML::convert<Property<std::string>>::decode(node[className], binaryFileBasename);
        unsignedProperty::decode(node[className], samplesToAnalyze);
        doubleProperty ::decode(node[className], minimalClusterWeight);
    }

    void MaximaProcessing::appendToNode(YAML::Node &node) const {
        node[className][binaryFileBasename.name()] = binaryFileBasename.get();
        node[className][samplesToAnalyze.name()] = samplesToAnalyze.get();
        node[className][minimalClusterWeight.name()] = minimalClusterWeight.get();
    }
}

YAML_SETTINGS_DEFINITION(Settings::MaximaProcessing)

namespace MaximaProcessing {
    Settings::MaximaProcessing settings = Settings::MaximaProcessing();
}