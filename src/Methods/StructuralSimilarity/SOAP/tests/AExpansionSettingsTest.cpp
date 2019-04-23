//
// Created by Michael Heuer on 2019-02-01.
//

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