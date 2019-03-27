//
// Created by Michael Heuer on 25.09.18.
//

#ifndef INPSIGHTS_GLOBALSIMILARITYSORTER_H
#define INPSIGHTS_GLOBALSIMILARITYSORTER_H

#include "SimilarReferences.h"
#include <HungarianHelper.h>
#include <ISettings.h>
#include <vector>

namespace Settings {
    class GlobalSimilaritySorter : public ISettings {
        inline static const std::string className = {VARNAME(GlobalSimilaritySorter)};
    public:
        Property<double> similarityRadius = {0.1, VARNAME(similarityRadius)};
        Property<double> similarityValueIncrement = {1e-5, VARNAME(similarityValueIncrement)};

        GlobalSimilaritySorter();
        explicit GlobalSimilaritySorter(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::GlobalSimilaritySorter)

class GlobalSimilaritySorter {
public:
    static Settings::GlobalSimilaritySorter settings;

    GlobalSimilaritySorter(
            std::vector<Sample> &samples,
            std::vector<Reference> &references,
            std::vector<SimilarReferences> &similarReferencesVector);
    bool sort();

private:
    std::vector<Sample>& samples_;
    std::vector<Reference>& references_;
    std::vector<SimilarReferences>& similarReferencesVector_;
};


#endif //INPSIGHTS_GLOBALSIMILARITYSORTER_H
