//
// Created by heuer on 12.12.18.
//

#ifndef INPSIGHTS_SOAPCLUSTERER_H
#define INPSIGHTS_SOAPCLUSTERER_H

#include "Sample.h"
#include "IClusterer.h"
#include <ISettings.h>
#include <MolecularSpectrum.h>

namespace Settings {
    class SOAPClusterer : public ISettings {
    public:
        Property<double> soapSimilarityThreshold = {0.95, VARNAME(soapSimilarityThreshold)};
        Property<double> distanceToleranceRadius = {0.2, VARNAME(distanceToleranceRadius)};
        SOAPClusterer();
        explicit SOAPClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::SOAPClusterer)

class SOAPClusterer : public IClusterer {
public:
    static Settings::SOAPClusterer settings;

    SOAPClusterer(const AtomsVector& atoms, std::vector<Sample> &samples);

    void cluster(Group &group) override;

private:
    AtomsVector atoms_;
    std::vector<Sample> &samples_;
};

#endif //INPSIGHTS_SOAPCLUSTERER_H
