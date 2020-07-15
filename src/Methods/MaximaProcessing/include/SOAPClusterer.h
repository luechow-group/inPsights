/* Copyright (C) 2018-2020 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_SOAPCLUSTERER_H
#define INPSIGHTS_SOAPCLUSTERER_H

#include "Sample.h"
#include "IBlock.h"
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
