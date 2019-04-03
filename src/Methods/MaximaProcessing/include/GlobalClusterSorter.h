//
// Created by heuer on 19.10.18.
//

#ifndef INPSIGHTS_GLOBALCLUSTERSORTER_H
#define INPSIGHTS_GLOBALCLUSTERSORTER_H

#include "SimilarReferences.h"
#include <DensityBasedScan.h>
#include <GlobalSimilaritySorter.h>
#include <ISettings.h>
#include <spdlog/spdlog.h>

namespace Settings {
    class GlobalClusterSorter : public ISettings {
        inline static const std::string className = {VARNAME(GlobalClusterSorter)};
    public:
        Property<double> clusterRadius = {
                ::GlobalSimilaritySorter::settings.similarityRadius(), VARNAME(clusterRadius)};

        GlobalClusterSorter();
        explicit GlobalClusterSorter(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::GlobalClusterSorter)

class GlobalClusterSorter {
public:
    static Settings::GlobalClusterSorter settings;

    GlobalClusterSorter(
            std::vector<Sample> &samples,
            std::vector<SimilarReferences> &globallySimilarMaxima,
            std::vector<std::vector<SimilarReferences>> &globallyClusteredMaxima);

    void sort();

private:
    std::vector<Sample> &samples_;
    std::vector<SimilarReferences> &globallySimilarMaxima_;
    std::vector<std::vector<SimilarReferences>> &globallyClusteredMaxima_;
    
    struct SortElement {
        SortElement(
                BestMatch::Result bestMatch,
                std::vector<SimilarReferences>::iterator it);

        bool operator<(const SortElement &rhs) const;

        BestMatch::Result bestMatch_;
        std::vector<SimilarReferences>::iterator it_;
    };

    static double wrapper(const SimilarReferences &s1, const SimilarReferences &s2);
};

#endif //INPSIGHTS_GLOBALCLUSTERSORTER_H
