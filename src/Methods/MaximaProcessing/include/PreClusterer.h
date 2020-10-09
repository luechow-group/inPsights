// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_PRECLUSTERER_H
#define INPSIGHTS_PRECLUSTERER_H

#include <BestMatch.h>
#include <IProcess.h>
#include <ISettings.h>

namespace Settings {
    class PreClusterer : public ISettings {
    public:
        Property<double> radius = {0.01, VARNAME(radius)};
        Property<double> valueIncrement = {1e-5, VARNAME(valueIncrement)};

        PreClusterer();

        explicit PreClusterer(const YAML::Node &node);

        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::PreClusterer)

class PreClusterer : public IClusterer {
public:
    static Settings::PreClusterer settings;

    PreClusterer(std::vector<Sample> &samples);

    void cluster(Cluster &group) override;
};


#endif //INPSIGHTS_PRECLUSTERER_H
