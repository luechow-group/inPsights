//
// Created by Michael Heuer on 28.08.2018.
//

#include <gmock/gmock.h>
#include "Hungarian.h"
#include "Metrics.h"
#include "TestMolecules.h"

TEST(HungarianTest, OneSolution) {
    Eigen::MatrixXd input(4,4);
    input << \
    88,83,69,92,\
    77,37,49,92,\
    11,69,5,86,\
    8,9,98,23;

    Eigen::VectorXi expectedOutput(4);
    expectedOutput << 2,1,0,3;

    ASSERT_EQ(Hungarian<double>::findMatching(input).indices(),expectedOutput);
}


TEST(HungarianTest, TwoSolutions) {
    Eigen::MatrixXd input(4,4);
    input <<\
    90,75,75,80,\
    35,85,55,65,\
    125,95,90,105,\
    45,110,95,115;

    Eigen::VectorXi expectedOutputs[2] = {Eigen::VectorXi(4),Eigen::VectorXi(4)};
    expectedOutputs[0]<< 1,3,2,0;
    expectedOutputs[1]<< 3,2,1,0;

    ASSERT_THAT(Hungarian<double>::findMatching(input).indices(),testing::AnyOf(expectedOutputs[0],expectedOutputs[1]));
}


TEST(HungarianTest, FlippedPositions) {
    auto p1 = TestMolecules::H2::ElectronsInCores::normal.electrons().positionsVector();;
    auto p2 = TestMolecules::H2::ElectronsInCores::flippedSpins.electrons().positionsVector();

    auto input = Metrics::positionalDistances(p1,p2);
    ASSERT_EQ(input.rows(),2);
    ASSERT_EQ(input.cols(),2);

    Eigen::VectorXi expectedOutput(2);
    expectedOutput << 1,0;
    ASSERT_EQ(Hungarian<double>::findMatching(input).indices(),expectedOutput);
}


TEST(HungarianTest, TranslatedPositions) {
    auto p1 = TestMolecules::H2::ElectronsInCores::normal.electrons().positionsVector();;
    auto p2 = TestMolecules::H2::ElectronsInCores::translated.electrons().positionsVector();

    auto input = Metrics::positionalDistances(p1,p2);
    ASSERT_EQ(input.rows(),2);
    ASSERT_EQ(input.cols(),2);

    Eigen::VectorXi expectedOutput(2);
    expectedOutput << 0,1;
    ASSERT_EQ(Hungarian<double>::findMatching(input).indices(),expectedOutput);
}

TEST(HungarianTest, TranslatedAndFlippedPositions) {

    auto p1 = TestMolecules::H2::ElectronsInCores::translated.electrons().positionsVector();;
    auto p2 = TestMolecules::H2::ElectronsInCores::flippedSpins.electrons().positionsVector();

    auto input = Metrics::positionalDistances(p1,p2);
    ASSERT_EQ(input.rows(),2);
    ASSERT_EQ(input.cols(),2);

    Eigen::VectorXi expectedOutput(2);
    expectedOutput << 1,0;
    ASSERT_EQ(Hungarian<double>::findMatching(input).indices(),expectedOutput);
}