/* Copyright (C) 2020 Michael Heuer.
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
#include <spdlog/spdlog.h>
#include <TotalWeightDifferenceAnalyzer.h>
#include <Reference.h>
#include <TestMolecules.h>

using namespace testing;

class ATotalWeightDifferenceAnalyzerTest : public ::testing::Test {
public:
    Group g1,g2,g3,g3a,g3b;

    void SetUp() override {
        g1 = Group({1,TestMolecules::H2::ElectronsInCores::normal.electrons(), 0});
        g2 = Group({1,TestMolecules::H2::ElectronsInCores::ionicRight.electrons(), 1});
        g3a = Group({2,TestMolecules::H2::ElectronsInCores::ionicLeft.electrons(), 2});
        g3b = Group({2,TestMolecules::H2::ElectronsInCores::ionicLeft.electrons(), 3});
        g3 = Group({g3a,g3b});
    }
};

TEST_F(ATotalWeightDifferenceAnalyzerTest, emptyGroupDeathTest){
    Group group({});

    TotalWeightDifferenceAnalyzer::settings.startRadius = 0.4;
    TotalWeightDifferenceAnalyzer::settings.increments = 3;
    TotalWeightDifferenceAnalyzer::settings.radiusIncrement = 0.2;  // Radii: before, 0.4, 0.6, 0.8, 1.2

    TotalWeightDifferenceAnalyzer analyzer;
    EXPECT_DEATH(analyzer.analyze(group),"");
}

TEST_F(ATotalWeightDifferenceAnalyzerTest, Ungrouped){
    Group group({g1,g2,g3a,g3b});


    std::vector<std::list<Eigen::Index>> previousClusterIndicesOfGroup;

    TotalWeightDifferenceAnalyzer::settings.startRadius = 0.4;
    TotalWeightDifferenceAnalyzer::settings.increments = 3;
    TotalWeightDifferenceAnalyzer::settings.radiusIncrement = 0.2; // Radii: before, 0.4, 0.6, 0.8, 1.2

    TotalWeightDifferenceAnalyzer analyzer;
    analyzer.analyze(group);

    ASSERT_THAT(analyzer.getResults(), ElementsAre(0.25, 0.0, 0.5, 0.0));
    //r=0.4: w(3a+3b)=0.5, w(3a)=0.25, Delta=0.25
    //r=0.6: w(1+2+3a+3b)=1.0, w(3a+3b)=0.5, Delta=0.5
}

TEST_F(ATotalWeightDifferenceAnalyzerTest, Pregrouped){
    Group group({g1,g2,g3});

    TotalWeightDifferenceAnalyzer::settings.startRadius = 0.4;
    TotalWeightDifferenceAnalyzer::settings.increments = 3;
    TotalWeightDifferenceAnalyzer::settings.radiusIncrement = 0.2; // Radii: before, 0.4, 0.6, 0.8, 1.2

    TotalWeightDifferenceAnalyzer analyzer;
    analyzer.analyze(group);

    ASSERT_THAT(analyzer.getResults(), ElementsAre(0.0, 0.0, 0.5, 0.0));
    //r=0.6: w(1+2+3a+3b)=1.0, w(3a+3b)=0.5, Delta=0.5
}
