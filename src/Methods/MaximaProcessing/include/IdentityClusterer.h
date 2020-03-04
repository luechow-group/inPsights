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

#ifndef INPSIGHTS_IDENTITYCLUSTERER_H
#define INPSIGHTS_IDENTITYCLUSTERER_H

#include "Sample.h"
#include "IBlock.h"
#include <ISettings.h>

namespace Settings {
    class IdentityClusterer : public ISettings {
    public:
        Property<double> radius = {0.01, VARNAME(radius)};
        Property<double> valueIncrement = {1e-7, VARNAME(valueIncrement)};

        IdentityClusterer();
        explicit IdentityClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::IdentityClusterer)

class IdentityClusterer : public IClusterer {
public:
    static Settings::IdentityClusterer settings;

    IdentityClusterer(std::vector<Sample> &samples);
    void cluster(Group& group) override;

private:
    void subLoop(Group& group,
            Group::iterator &beginIt,
            Group::iterator &it,
            Group::iterator &endIt,
            double distThresh,
            double valueIncrement);

    // TODO This method should be located inside of a reference container class
    void addReference(Group& group,
            const Group::iterator &beginIt,
            Group::iterator &it,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const;
};


#endif //INPSIGHTS_IDENTITYCLUSTERER_H
