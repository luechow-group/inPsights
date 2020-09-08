// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
    void cluster(Cluster& group) override;

private:
    void subLoop(Cluster& group,
            Cluster::iterator &beginIt,
            Cluster::iterator &it,
            Cluster::iterator &endIt,
            double distThresh,
            double valueIncrement);

    // TODO This method should be located inside of a reference container class
    void addReference(Cluster& group,
            const Cluster::iterator &beginIt,
            Cluster::iterator &it,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const;
};


#endif //INPSIGHTS_IDENTITYCLUSTERER_H
