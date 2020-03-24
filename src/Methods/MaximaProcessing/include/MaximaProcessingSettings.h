/* Copyright (C) 2018-2019 Michael Heuer.
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

#ifndef INPSIGHTS_MAXIMAPROCESSINGSETTINGS_H
#define INPSIGHTS_MAXIMAPROCESSINGSETTINGS_H

#include <ISettings.h>
#include <Property.h>
#include <string>

enum class SamplesToAnalyze{
    small = 10000,
    medium = 100000,
    large = 1000000,
    all = 0 };

namespace Settings {
    class MaximaProcessing : public ISettings {
    public:
        Property<std::string> binaryFileBasename = {"", VARNAME(binaryFileBasename)};
        Property<unsigned> samplesToAnalyze = {unsigned(SamplesToAnalyze::small), VARNAME(samplesToAnalyze)};
        Property<double> minimalClusterWeight = {0.01, VARNAME(minimalClusterWeight)};
        Property<bool> deleteCoreElectrons = {false, VARNAME(deleteCoreElectrons)};

        MaximaProcessing();
        explicit MaximaProcessing(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::MaximaProcessing)

namespace MaximaProcessing {
    extern Settings::MaximaProcessing settings;
}

#endif //INPSIGHTS_MAXIMAPROCESSINGSETTINGS_H
