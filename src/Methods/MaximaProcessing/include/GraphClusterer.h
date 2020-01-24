/* Copyright 2020 heuer
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

#ifndef INPSIGHTS_GRAPHCLUSTERER_H
#define INPSIGHTS_GRAPHCLUSTERER_H

#include "Sample.h"
#include "IClusterer.h"
#include <ISettings.h>

namespace Settings {
    class GraphClusterer : public ISettings {
    public:
        Property<double> startRadius = {0.0, VARNAME(startRadius)};
        Property<double> endRadius = {1.0, VARNAME(endRadius)};
        Property<double> radiusIncrement = {0.05, VARNAME(radiusIncrement)};

        GraphClusterer();
        explicit GraphClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::GraphClusterer)
class GraphClusterer : public IClusterer{
public:
    static Settings::GraphClusterer settings;

    GraphClusterer(Group& group);
    Eigen::MatrixXd calculateAdjacencyMatrix(Group& group);
    void cluster(Group& group) override;
    std::vector<std::size_t> scanClusterSizeWithDistance();

private:
    //std::vector<Sample> &samples_;
    Eigen::MatrixXd mat_;
};

#endif //INPSIGHTS_GRAPHCLUSTERER_H
