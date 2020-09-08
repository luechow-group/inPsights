// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_TOTALWEIGHTDIFFERENCEANALYZER_H
#define INPSIGHTS_TOTALWEIGHTDIFFERENCEANALYZER_H

#include "IBlock.h"
#include <ISettings.h>

namespace Settings {
    class TotalWeightDifferenceAnalyzer : public ISettings {
    public:
        Property<double> startRadius = {0.0, VARNAME(startRadius)};
        Property<unsigned> increments = {30, VARNAME(increments)};
        Property<double> radiusIncrement = {0.05, VARNAME(radiusIncrement)};

        TotalWeightDifferenceAnalyzer();
        explicit TotalWeightDifferenceAnalyzer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::TotalWeightDifferenceAnalyzer)
class TotalWeightDifferenceAnalyzer : public IAnalyzer {
public:
    static Settings::TotalWeightDifferenceAnalyzer settings;

    void analyze(const Cluster& group) override;
    std::vector<double> getResults();

private:
    std::vector<double> totalWeightDifferences_;
};


#endif //INPSIGHTS_TOTALWEIGHTDIFFERENCEANALYZER_H
