// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_DENSITYBASEDCLUSTERER_H
#define INPSIGHTS_DENSITYBASEDCLUSTERER_H

#include <DensityBasedScan.h>
#include <PreClusterer.h>
#include <ISettings.h>
#include <spdlog/spdlog.h>
#include <IBlock.h>
#include <Cluster.h>

namespace Settings {
    class DensityBasedClusterer : public ISettings {
    public:
        Property<double> radius = {0.2, VARNAME(radius)};
        Property<size_t> minimalClusterSize = {1, VARNAME(minimalClusterSize)};
        Property<bool> local =  {false, VARNAME(local)}; // TODO unite all general clusterer settigns in a parent settings class
        Property<bool> sortRemainder = {false, VARNAME(sortRemainder)};

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

    void cluster(Cluster& group) override;

private:
    static double wrapper(const Cluster &g1, const Cluster &g2);
    static double wrapperLocal(const Cluster &g1, const Cluster &g2);
    void orderByBestMatchDistance(Cluster &supergroup, double threshold, bool local = false) const;

    bool compare(double threshold, Cluster &subgroup, Cluster &newClusters,
            const std::vector<Cluster>::iterator &i,
            std::vector<Cluster>::iterator &j) const;

    bool compareLocal(double threshold, Cluster &subgroup, Cluster &newClusters,
            const std::vector<Cluster>::iterator &i,
            std::vector<Cluster>::iterator &j) const;
};

#endif //INPSIGHTS_DENSITYBASEDCLUSTERER_H
