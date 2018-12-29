//
// Created by Michael Heuer on 2018-12-28.
//

#include "GlobalSortSettings.h"
#include <assert.h>
#include <Logger.h>

namespace Settings {
    GlobalSort::GlobalSort() {
        samplesToAnalyze.onChange().connect([&](unsigned value) {
            if(value > 0)
                Logger::console->info("Analyzing {} samples.", value);
            else if (value == 0) {
                samplesToAnalyze = std::numeric_limits<unsigned>::max();
                Logger::console->info("Analyzing all samples.");
            } else
                throw std::invalid_argument("The number of samples to analyze is negative.");
        });

        identityRadius.onChange().connect([&](double value) {
            if(value > similarityRadius.get())
                throw std::invalid_argument("Choosing the identityRadius greater than the similarityRadius makes no sense.");
        });
        similarityRadius.onChange().connect([&](double value) {
            if(value < identityRadius.get())
                throw std::invalid_argument("Choosing the similarityRadius smaller than the identityRadius makes no sense.");

        });
    }

    GlobalSort::GlobalSort(const YAML::Node &node)
    : GlobalSort() {
        unsignedProperty::decode(node[className], samplesToAnalyze);
        boolProperty::decode(node[className], identitySearch);
        doubleProperty::decode(node[className], identityRadius);
        doubleProperty::decode(node[className], similarityRadius);
    }

    void GlobalSort::addToNode(YAML::Node &node) const {
        node[className][samplesToAnalyze.name()] = samplesToAnalyze.get();
        node[className][identitySearch.name()] = identitySearch.get();
        node[className][identityRadius.name()] = identityRadius.get();
        node[className][similarityRadius.name()] = similarityRadius.get();
    }
}

YAML_SETTINGS_DEFINITION(Settings::GlobalSort)