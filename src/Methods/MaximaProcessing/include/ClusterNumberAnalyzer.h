/* Copyright (C) 2020 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
