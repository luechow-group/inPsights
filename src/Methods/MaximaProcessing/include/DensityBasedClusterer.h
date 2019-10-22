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

#ifndef INPSIGHTS_DENSITYBASEDCLUSTERER_H
#define INPSIGHTS_DENSITYBASEDCLUSTERER_H

#include <DensityBasedScan.h>
#include <DistanceClusterer.h>
#include <ISettings.h>
#include <spdlog/spdlog.h>

#include <IClusterer.h>
#include <Group.h>

namespace Settings {
    class DensityBasedClusterer : public ISettings {
    public:
        Property<double> radius = {::DistanceClusterer::settings.radius(), VARNAME(radius)};
        Property<size_t> minimalClusterSize = {1, VARNAME(minimalClusterSize)};

        DensityBasedClusterer();
        explicit DensityBasedClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::DensityBasedClusterer)


class DensityBasedClusterer : public IClusterer {
public:
    static Settings::DensityBasedClusterer settings;

    explicit DensityBasedClusterer(std::vector<Sample> &samples);

    void cluster(Group& group) override;

private:
    std::vector<Sample> &samples_;

    static double wrapper(const Group &g1, const Group &g2);

    void orderByBestMatchDistance(Group &supergroup, double threshold) const;
};

#endif //INPSIGHTS_DENSITYBASEDCLUSTERER_H
