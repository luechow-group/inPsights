/* Copyright (C) 2019 Leonard Reuter.
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

#ifndef INPSIGHTS_REFERENCEPOSITIONSCLUSTERER_H
#define INPSIGHTS_REFERENCEPOSITIONSCLUSTERER_H

#include <IBlock.h>
#include <ISettings.h>
#include <ParticlesVector.h>
#include <Sample.h>
#include <functional>

namespace Settings {
    class ReferencePositionsClusterer : public ISettings {
    public:
        Property<double> radius = {0.1, VARNAME(radius)};
        Property<bool> local =  {false, VARNAME(local)}; // TODO unite all general clusterer settigns in a parent settings class
        Property<bool> sortRemainder = {false, VARNAME(sortRemainder)};

        ReferencePositionsClusterer();
        explicit ReferencePositionsClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::ReferencePositionsClusterer)


class ReferencePositionsClusterer : public IClusterer{
public:
    static Settings::ReferencePositionsClusterer settings;

    ReferencePositionsClusterer(std::vector<Sample> &samples);

    void cluster(Cluster& group) override;

private:
    std::function<double(const Eigen::Vector3d &, const std::vector<Eigen::Vector3d> &)> distanceFunction_;
    bool compareLocal(std::vector<Cluster>::iterator &sortedGroup, std::vector<Cluster>::iterator &subGroup,
            double similarityRadius) const;

    bool compareGlobal(std::vector<Cluster>::iterator &sortedGroup, std::vector<Cluster>::iterator &subGroup,
            double similarityRadius) const;
};

#endif //INPSIGHTS_REFERENCEPOSITIONSCLUSTERER_H
