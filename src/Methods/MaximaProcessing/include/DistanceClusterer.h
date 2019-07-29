//
// Created by Michael Heuer on 25.09.18.
//

#ifndef INPSIGHTS_DISTANCECLUSTERER_H
#define INPSIGHTS_DISTANCECLUSTERER_H

#include <BestMatch.h>
#include <IClusterer.h>
#include <ISettings.h>

namespace Settings {
    class DistanceClusterer : public ISettings {
    public:
        Property<double> radius = {0.1, VARNAME(radius)};
        Property<double> valueIncrement = {1e-5, VARNAME(valueIncrement)};

        DistanceClusterer();
        explicit DistanceClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::DistanceClusterer)

class DistanceClusterer : public IClusterer {
public:
    static Settings::DistanceClusterer settings;

    DistanceClusterer(std::vector<Sample> &samples);
    void cluster(Group& group) override;

private:
    std::vector<Sample>& samples_;
};


#endif //INPSIGHTS_DISTANCECLUSTERER_H
