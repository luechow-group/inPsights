//
// Created by Michael Heuer on 25.09.18.
//

#ifndef INPSIGHTS_BESTMATCHDISTANCESIMILARITYCLUSTERER_H
#define INPSIGHTS_BESTMATCHDISTANCESIMILARITYCLUSTERER_H

#include <BestMatch.h>
#include <IClusterer.h>
#include <ISettings.h>

namespace Settings {
    class BestMatchDistanceSimilarityClusterer : public ISettings {
    public:
        Property<double> similarityRadius = {0.1, VARNAME(similarityRadius)};
        Property<double> similarityValueIncrement = {1e-5, VARNAME(similarityValueIncrement)};

        BestMatchDistanceSimilarityClusterer();
        explicit BestMatchDistanceSimilarityClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::BestMatchDistanceSimilarityClusterer)

class BestMatchDistanceSimilarityClusterer : public IClusterer {
public:
    static Settings::BestMatchDistanceSimilarityClusterer settings;

    BestMatchDistanceSimilarityClusterer(std::vector<Sample> &samples);
    void cluster(Group& group) override;

private:
    std::vector<Sample>& samples_;
};


#endif //INPSIGHTS_BESTMATCHDISTANCESIMILARITYCLUSTERER_H
