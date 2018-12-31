//
// Created by Michael Heuer on 2018-12-28.
//

#ifndef INPSIGHTS_ENERGYPARTITIONINGSETTINGS_H
#define INPSIGHTS_ENERGYPARTITIONINGSETTINGS_H

#include <ISettings.h>
#include <Property.h>
#include <string>

enum class SamplesToAnalyze{
    small = 10000,
    medium = 100000,
    large = 1000000,
    all = 0 };

namespace Settings {
    class EnergyPartitioning : public ISettings {
    public:
        inline static const std::string className = {VARNAME(EnergyPartitioning)};

        Property<std::string> binaryFileBasename = {"raw", VARNAME(binaryFileBasename)};
        Property<unsigned> samplesToAnalyze = {unsigned(SamplesToAnalyze::small), VARNAME(samplesToAnalyze)};
        Property<bool> identitySearch = {false, VARNAME(identitySearch)};

        EnergyPartitioning();
        explicit EnergyPartitioning(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::EnergyPartitioning)

#endif //INPSIGHTS_ENERGYPARTITIONINGSETTINGS_H
