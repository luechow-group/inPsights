//
// Created by Michael Heuer on 2018-12-28.
//

#include "GlobalSortSettings.h"
#include <assert.h>

namespace Settings {
    GlobalSort::GlobalSort() {
        identityRadius.onChange().connect([&](double val) {
            assert(val > 0 && "The threshold cannot be negative.");
        });
    }

    GlobalSort::GlobalSort(const YAML::Node &node) {
        intProperty::decode(node[className], samplesToAnalyze);
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