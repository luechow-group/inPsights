//
// Created by heuer on 19.10.18.
//

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
        Property<double> radius = {
                ::DistanceClusterer::settings.radius(), VARNAME(radius)};

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
