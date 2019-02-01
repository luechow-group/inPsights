//
// Created by Michael Heuer on 2019-02-01.
//

#include <gmock/gmock.h>
#include <ExpansionSettings.h>
#include <sstream>

TEST(AExpansionSettingsTest, YAMLConversion){
    YAML::Node node;
    SOAPExpansion::settings.appendToNode(node);
    Radial::settings.appendToNode(node);
    Angular::settings.appendToNode(node);
    Cutoff::settings.appendToNode(node);

    std::stringstream ss;
    ss << node;

    auto decodedSettings = Settings::SOAPExpansion(node);

    auto decodedNode = YAML::Load(ss.str().c_str());
    auto decodedSettings2 = Settings::SOAPExpansion(decodedNode);

    auto refPairSim = SOAPExpansion::settings.pairSimilarities[{-2,-1}];
    auto pairSim = decodedSettings.pairSimilarities[{-2,-1}];
    auto pairSim2 = decodedSettings.pairSimilarities[{-2,-1}];
    ASSERT_EQ(pairSim, refPairSim);
    ASSERT_EQ(pairSim2, refPairSim);
}