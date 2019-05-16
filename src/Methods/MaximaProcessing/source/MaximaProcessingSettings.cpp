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
                    else if (value < minimalClusterSize.get())
                        throw std::invalid_argument("The number of samples cannot be smaller than the minimal cluster size");
                });
        minimalClusterSize.onChange_.connect(
                [&](unsigned value) {
                    if (value > samplesToAnalyze.get())
                        throw std::invalid_argument("The minimal cluster size cannot be greater than the number of samples");
                });
    }

    MaximaProcessing::MaximaProcessing(const YAML::Node &node)
    : MaximaProcessing() {
        YAML::convert<Property<std::string>>::decode(node[className], binaryFileBasename);
        unsignedProperty::decode(node[className], samplesToAnalyze);
        unsignedProperty::decode(node[className], minimalClusterSize);
    }

    void MaximaProcessing::appendToNode(YAML::Node &node) const {
        node[className][binaryFileBasename.name()] = binaryFileBasename.get();
        node[className][samplesToAnalyze.name()] = samplesToAnalyze.get();
        node[className][minimalClusterSize.name()] = minimalClusterSize.get();
    }
}

YAML_SETTINGS_DEFINITION(Settings::MaximaProcessing)

namespace MaximaProcessing {
    Settings::MaximaProcessing settings = Settings::MaximaProcessing();
}