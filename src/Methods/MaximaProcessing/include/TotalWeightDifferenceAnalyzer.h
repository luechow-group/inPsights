/* Copyright 2020 heuer
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

#ifndef INPSIGHTS_TOTALWEIGHTDIFFERENCEANALYZER_H
#define INPSIGHTS_TOTALWEIGHTDIFFERENCEANALYZER_H

#include "IBlock.h"
#include <ISettings.h>

namespace Settings {
    class TotalWeightDifferenceAnalyzer : public ISettings {
    public:
        Property<double> startRadius = {0.0, VARNAME(startRadius)};
        Property<unsigned> increments = {30, VARNAME(increments)};
        Property<double> radiusIncrement = {0.05, VARNAME(radiusIncrement)};

        TotalWeightDifferenceAnalyzer();
        explicit TotalWeightDifferenceAnalyzer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::TotalWeightDifferenceAnalyzer)
class TotalWeightDifferenceAnalyzer : public IAnalyzer {
public:
    static Settings::TotalWeightDifferenceAnalyzer settings;

    void analyze(const Group& group) override;
    std::vector<double> getResults();

private:
    std::vector<double> totalWeightDifferences_;
};

template<typename K, typename V>
bool findByValue(std::vector<K> & vec, std::map<K, V> mapOfElement, V value)
{
    bool bResult = false;
    auto it = mapOfElement.begin();
    // Iterate through the map
    while(it != mapOfElement.end())
    {
        // Check if value of this entry matches with given value
        if(it->second == value)
        {
            // Yes found
            bResult = true;
            // Push the key in given map
            vec.push_back(it->first);
        }
        // Go to next entry in map
        it++;
    }
    return bResult;
}

#endif //INPSIGHTS_TOTALWEIGHTDIFFERENCEANALYZER_H
