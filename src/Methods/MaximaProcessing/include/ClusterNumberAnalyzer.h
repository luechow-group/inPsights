// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CLUSTERNUMBERANALYZER_H
#define INPSIGHTS_CLUSTERNUMBERANALYZER_H

#include "IBlock.h"
#include <ISettings.h>

namespace Settings {
    class ClusterNumberAnalyzer : public ISettings {
    public:
        Property<double> startRadius = {0.0, VARNAME(startRadius)};
        Property<double> radiusIncrement = {0.05, VARNAME(radiusIncrement)};
        Property<unsigned> increments = {30, VARNAME(increments)};
        Property<double> minimalWeight = {0.0, VARNAME(minimalWeight)};

        ClusterNumberAnalyzer();
        explicit ClusterNumberAnalyzer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::ClusterNumberAnalyzer)
class ClusterNumberAnalyzer : public IAnalyzer {
public:
    static Settings::ClusterNumberAnalyzer settings;

    void analyze(const Cluster& cluster) override;
    std::vector<std::size_t> getResults();

private:
    std::vector<std::size_t> clusterNumbers_;
};


#endif //INPSIGHTS_CLUSTERNUMBERANALYZER_H
