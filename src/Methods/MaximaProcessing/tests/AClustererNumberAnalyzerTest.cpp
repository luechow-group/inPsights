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
#include <ClusterNumberAnalyzer.h>
#include <Reference.h>
#include <TestMolecules.h>

using namespace testing;

class AClusterNumberAnalyzerTest : public ::testing::Test {
public:
    Cluster g1,g2,g3,g3a,g3b;
    AtomsVector atoms;

    void SetUp() override {
        atoms = TestMolecules::H2::ElectronsInCores::normal.atoms();

        g1 = Cluster({atoms, 1,TestMolecules::H2::ElectronsInCores::normal.electrons(), 0});
        g2 = Cluster({atoms, 1,TestMolecules::H2::ElectronsInCores::ionicRight.electrons(), 1});
        g3a = Cluster({atoms, 2,TestMolecules::H2::ElectronsInCores::ionicLeft.electrons(), 2});
        g3b = Cluster({atoms, 2,TestMolecules::H2::ElectronsInCores::ionicLeft.electrons(), 3});
        g3 = Cluster({g3a,g3b});
    }
};

TEST_F(AClusterNumberAnalyzerTest, DISABLED_emptyClusterDeathTest){
    Cluster group;

    ClusterNumberAnalyzer::settings.startRadius = 0.4;
    ClusterNumberAnalyzer::settings.increments = 2;
    ClusterNumberAnalyzer::settings.radiusIncrement = 0.2;  // Radii: before, 0.4, 0.6, 0.8

    ClusterNumberAnalyzer analyzer;
    EXPECT_DEATH(analyzer.analyze(group),"");
}

TEST_F(AClusterNumberAnalyzerTest, normalUse){
    Cluster group({g1,g2,g3});

    ClusterNumberAnalyzer::settings.startRadius = 0.4;
    ClusterNumberAnalyzer::settings.increments = 2;
    ClusterNumberAnalyzer::settings.radiusIncrement = 0.2; // Radii: before, 0.4, 0.6, 0.8

    ClusterNumberAnalyzer analyzer;
    analyzer.analyze(group);

    ASSERT_THAT(analyzer.getResults(), ElementsAre(3,3,3,1));
}

TEST_F(AClusterNumberAnalyzerTest, minimalWeight){
    Cluster group({g1,g2,g3});

    ClusterNumberAnalyzer::settings.startRadius = 0.4;
    ClusterNumberAnalyzer::settings.increments = 2;
    ClusterNumberAnalyzer::settings.radiusIncrement = 0.2;
    ClusterNumberAnalyzer::settings.minimalWeight = 0.3;  // Radii: before, 0.4, 0.6, 0.8

    ClusterNumberAnalyzer analyzer;
    analyzer.analyze(group);

    ASSERT_THAT(analyzer.getResults(), ElementsAre(1,1,1,1));
}