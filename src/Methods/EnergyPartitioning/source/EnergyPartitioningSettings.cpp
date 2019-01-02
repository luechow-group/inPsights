//
// Created by Michael Heuer on 2018-12-28.
//

#include "EnergyPartitioningSettings.h"
#include <assert.h>
#include <Logger.h>

namespace Settings {
    EnergyPartitioning::EnergyPartitioning() {
        samplesToAnalyze.onChange().connect([&](unsigned value) {
            if(value > 0 && value < std::numeric_limits<unsigned>::max())
                Logger::console->info("Analyzing {} samples.", value);
            else if (value == 0)
                samplesToAnalyze = std::numeric_limits<unsigned>::max();
            else if (value == std::numeric_limits<unsigned>::max())
                Logger::console->info("Analyzing all samples.");
            else
                throw std::invalid_argument("The number of samples to analyze is negative.");
        });
    }

    EnergyPartitioning::EnergyPartitioning(const YAML::Node &node)
    : EnergyPartitioning() {
        YAML::convert<Property<std::string>>::decode(node[className], binaryFileBasename);
        unsignedProperty::decode(node[className], samplesToAnalyze);
        boolProperty::decode(node[className], identitySearch);
    }

    void EnergyPartitioning::appendToNode(YAML::Node &node) const {
        node[className][binaryFileBasename.name()] = binaryFileBasename.get();
        node[className][samplesToAnalyze.name()] = samplesToAnalyze.get();
        node[className][identitySearch.name()] = identitySearch.get();
    }
}
YAML_SETTINGS_DEFINITION(Settings::EnergyPartitioning)
