//
// Created by Michael Heuer on 2018-12-28.
//

#include "MaximaProcessingSettings.h"
#include <assert.h>
#include <spdlog/spdlog.h>

namespace Settings {
    MaximaProcessing::MaximaProcessing() {
        samplesToAnalyze.onChange_.connect(
                [&](unsigned value) {
                    if(value > 0 && value < std::numeric_limits<unsigned>::max())
                        spdlog::info("Analyzing {} samples.", value);
                    else if (value == 0)
                        spdlog::info("Analyzing all samples.");
                    else
                        throw std::invalid_argument("The number of samples to analyze is negative.");
                });
    }


    MaximaProcessing::MaximaProcessing(const YAML::Node &node)
    : MaximaProcessing() {
        YAML::convert<Property<std::string>>::decode(node[className], binaryFileBasename);
        unsignedProperty::decode(node[className], samplesToAnalyze);
        boolProperty::decode(node[className], identitySearch);
    }

    void MaximaProcessing::appendToNode(YAML::Node &node) const {
        node[className][binaryFileBasename.name()] = binaryFileBasename.get();
        node[className][samplesToAnalyze.name()] = samplesToAnalyze.get();
        node[className][identitySearch.name()] = identitySearch.get();
    }
}
YAML_SETTINGS_DEFINITION(Settings::MaximaProcessing)
