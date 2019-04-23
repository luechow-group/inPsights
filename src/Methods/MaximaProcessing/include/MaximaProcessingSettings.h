//
// Created by Michael Heuer on 2018-12-28.
//

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
        Property<std::string> binaryFileBasename = {"raw", VARNAME(binaryFileBasename)};
        Property<unsigned> samplesToAnalyze = {unsigned(SamplesToAnalyze::small), VARNAME(samplesToAnalyze)};

        MaximaProcessing();
        explicit MaximaProcessing(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::MaximaProcessing)

#endif //INPSIGHTS_MAXIMAPROCESSINGSETTINGS_H
