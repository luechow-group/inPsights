/* Copyright (C) 2019 Michael Heuer.
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

#include <gmock/gmock.h>
#include <SOAPSettings.h>
#include <sstream>

using namespace SOAP;

TEST(AExpansionSettingsTest, YAMLConversion){
    YAML::Node node;
    General::settings.appendToNode(node);
    Radial::settings.appendToNode(node);
    Angular::settings.appendToNode(node);
    Cutoff::settings.appendToNode(node);

    std::stringstream ss;
    ss << node;

    auto decodedSettings = Settings::SOAP::General(node);

    auto decodedNode = YAML::Load(ss.str().c_str());
    auto decodedSettings2 = Settings::SOAP::General(decodedNode);

    auto refPairSim = General::settings.pairSimilarities[{-2,-1}];
    auto pairSim = decodedSettings.pairSimilarities[{-2,-1}];
    auto pairSim2 = decodedSettings.pairSimilarities[{-2,-1}];
    ASSERT_EQ(pairSim, refPairSim);
    ASSERT_EQ(pairSim2, refPairSim);
}