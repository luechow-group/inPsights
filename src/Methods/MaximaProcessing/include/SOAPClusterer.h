// Copyright (C) 2018-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SOAPCLUSTERER_H
#define INPSIGHTS_SOAPCLUSTERER_H

#include "Sample.h"
#include "IProcess.h"
#include <ISettings.h>
#include <MolecularSpectrum.h>

namespace Settings {
    class SOAPClusterer : public ISettings {
    public:
        Property<double> similarityThreshold = {0.98, VARNAME(similarityThreshold)};
        Property<double> distanceMatrixCovarianceTolerance = {0.2, VARNAME(distanceMatrixCovarianceTolerance)};
        Property<double> maxValueDelta = {1e-2, VARNAME(maxValueDelta)};
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

    void cluster(Cluster &cluster) override;

private:
    AtomsVector atoms_;
};

#endif //INPSIGHTS_SOAPCLUSTERER_H
