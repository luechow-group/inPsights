// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <spdlog/spdlog.h>
#include <TotalWeightDifferenceAnalyzer.h>
#include <Reference.h>
#include <TestMolecules.h>

using namespace testing;

class ATotalWeightDifferenceAnalyzerTest : public ::testing::Test {
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

TEST_F(ATotalWeightDifferenceAnalyzerTest, DISABLED_emptyGroupDeathTest){
    Cluster cluster;

    TotalWeightDifferenceAnalyzer::settings.startRadius = 0.4;
    TotalWeightDifferenceAnalyzer::settings.increments = 3;
    TotalWeightDifferenceAnalyzer::settings.radiusIncrement = 0.2;  // Radii: before, 0.4, 0.6, 0.8, 1.2

    TotalWeightDifferenceAnalyzer analyzer;
    EXPECT_DEATH(analyzer.analyze(cluster),"");
}

TEST_F(ATotalWeightDifferenceAnalyzerTest, Unclustered){
    Cluster cluster({g1,g2,g3a,g3b});


    std::vector<std::list<Eigen::Index>> previousClusterIndicesOfGroup;

    TotalWeightDifferenceAnalyzer::settings.startRadius = 0.4;
    TotalWeightDifferenceAnalyzer::settings.increments = 3;
    TotalWeightDifferenceAnalyzer::settings.radiusIncrement = 0.2; // Radii: before, 0.4, 0.6, 0.8, 1.2

    TotalWeightDifferenceAnalyzer analyzer;
    analyzer.analyze(cluster);

    ASSERT_THAT(analyzer.getResults(), ElementsAre(0.25, 0.0, 0.5, 0.0));
    //r=0.4: w(3a+3b)=0.5, w(3a)=0.25, Delta=0.25
    //r=0.6: w(1+2+3a+3b)=1.0, w(3a+3b)=0.5, Delta=0.5
}

TEST_F(ATotalWeightDifferenceAnalyzerTest, Preclustered){
    Cluster cluster({g1,g2,g3});

    TotalWeightDifferenceAnalyzer::settings.startRadius = 0.4;
    TotalWeightDifferenceAnalyzer::settings.increments = 3;
    TotalWeightDifferenceAnalyzer::settings.radiusIncrement = 0.2; // Radii: before, 0.4, 0.6, 0.8, 1.2

    TotalWeightDifferenceAnalyzer analyzer;
    analyzer.analyze(cluster);

    ASSERT_THAT(analyzer.getResults(), ElementsAre(0.0, 0.0, 0.5, 0.0));
    //r=0.6: w(1+2+3a+3b)=1.0, w(3a+3b)=0.5, Delta=0.5
}
