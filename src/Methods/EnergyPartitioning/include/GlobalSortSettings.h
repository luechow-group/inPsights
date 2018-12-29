//
// Created by Michael Heuer on 2018-12-28.
//

#ifndef INPSIGHTS_GLOBALSORTSETTINGS_H
#define INPSIGHTS_GLOBALSORTSETTINGS_H

#include <ISettings.h>
#include <Property.h>


enum class SamplesToAnalyze{
    small = 10000,
    medium = 100000,
    large = 1000000,
    all = -1 };

namespace Settings {
    class GlobalSort : public ISettings {
    public:
        inline static const std::string className = {VARNAME(GlobalSort)};
        Property<int> samplesToAnalyze = {int(SamplesToAnalyze::small), VARNAME(samplesToAnalyze)};
        Property<bool> identitySearch = {false, VARNAME(identitySearch)};
        Property<double> identityRadius = {0.01, VARNAME(identityRadius)};
        Property<double> similarityRadius = {0.1, VARNAME(similarityRadius)};

        GlobalSort();
        explicit GlobalSort(const YAML::Node &node);
        void addToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::GlobalSort)


#endif //INPSIGHTS_GLOBALSORTSETTINGS_H
