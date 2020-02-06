/* Copyright (C) 2018-2019 Michael Heuer.
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

#ifndef INPSIGHTS_DISTANCECLUSTERER_H
#define INPSIGHTS_DISTANCECLUSTERER_H

#include <BestMatch.h>
#include <IBlock.h>
#include <ISettings.h>

namespace Settings {
    class DistanceClusterer : public ISettings {
    public:
        Property<double> radius = {0.01, VARNAME(radius)};
        Property<double> valueIncrement = {1e-5, VARNAME(valueIncrement)};

        DistanceClusterer();
        explicit DistanceClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::DistanceClusterer)

class DistanceClusterer : public IClusterer {
public:
    static Settings::DistanceClusterer settings;

    DistanceClusterer(std::vector<Sample> &samples);
    void cluster(Group& group) override;

private:
    std::vector<Sample>& samples_;
};


#endif //INPSIGHTS_DISTANCECLUSTERER_H
