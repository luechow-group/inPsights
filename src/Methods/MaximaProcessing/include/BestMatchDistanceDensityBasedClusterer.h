//
// Created by heuer on 19.10.18.
//

#ifndef INPSIGHTS_BESTMATCHDISTANCEDENSITYBASEDCLUSTERER_H
#define INPSIGHTS_BESTMATCHDISTANCEDENSITYBASEDCLUSTERER_H

#include <DensityBasedScan.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <ISettings.h>
#include <spdlog/spdlog.h>

#include <IClusterer.h>
#include <Group.h>

namespace Settings {
    class BestMatchDistanceDensityBasedClusterer : public ISettings { //TODO rename
    public:
        Property<double> clusterRadius = {
                ::BestMatchDistanceSimilarityClusterer::settings.similarityRadius(), VARNAME(clusterRadius)};

        BestMatchDistanceDensityBasedClusterer();
        explicit BestMatchDistanceDensityBasedClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::BestMatchDistanceDensityBasedClusterer)


class BestMatchDistanceDensityBasedClusterer : public IClusterer {
public:
    static Settings::BestMatchDistanceDensityBasedClusterer settings;

    explicit BestMatchDistanceDensityBasedClusterer(std::vector<Sample> &samples);

    void cluster(Group& group) override;

private:
    std::vector<Sample> &samples_;

    static double wrapper(const Group &g1, const Group &g2);

    struct SortElement {
        SortElement(BestMatch::Result bestMatch, Group::iterator it);

        bool operator<(const SortElement &rhs) const;

        BestMatch::Result bestMatch_;
        Group::iterator it_;
    };

    void orderByBestMatchDistance(Group &supergroup, double threshold) const;
};

#endif //INPSIGHTS_BESTMATCHDISTANCEDENSITYBASEDCLUSTERER_H
