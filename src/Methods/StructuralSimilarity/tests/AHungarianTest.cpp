//
// Created by Michael Heuer on 28.08.2018.
//

#include <gmock/gmock.h>
#include <Hungarian.h>
#include <Metrics.h>
#include <TestMolecules.h>
#include <Interval.h>
#include <algorithm>
#include <random>

TEST(AHungarianTest, OneSolution) {
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


TEST(AHungarianTest, TwoSolutions) {
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


TEST(AHungarianTest, FlippedPositions) {
    auto p1 = TestMolecules::H2::ElectronsInCores::normal.electrons().positionsVector();;
    auto p2 = TestMolecules::H2::ElectronsInCores::flippedSpins.electrons().positionsVector();

    auto input = Metrics::positionalDistances(p1,p2);
    ASSERT_EQ(input.rows(),2);
    ASSERT_EQ(input.cols(),2);

    Eigen::VectorXi expectedOutput(2);
    expectedOutput << 1,0;
    ASSERT_EQ(Hungarian<double>::findMatching(input).indices(),expectedOutput);
}


TEST(AHungarianTest, TranslatedPositions) {
    auto p1 = TestMolecules::H2::ElectronsInCores::normal.electrons().positionsVector();;
    auto p2 = TestMolecules::H2::ElectronsInCores::translated.electrons().positionsVector();

    auto input = Metrics::positionalDistances(p1,p2);
    ASSERT_EQ(input.rows(),2);
    ASSERT_EQ(input.cols(),2);

    Eigen::VectorXi expectedOutput(2);
    expectedOutput << 0,1;
    ASSERT_EQ(Hungarian<double>::findMatching(input).indices(),expectedOutput);
}

TEST(AHungarianTest, TranslatedAndFlippedPositions) {

    auto p1 = TestMolecules::H2::ElectronsInCores::translated.electrons().positionsVector();;
    auto p2 = TestMolecules::H2::ElectronsInCores::flippedSpins.electrons().positionsVector();

    auto input = Metrics::positionalDistances(p1,p2);
    ASSERT_EQ(input.rows(),2);
    ASSERT_EQ(input.cols(),2);

    Eigen::VectorXi expectedOutput(2);
    expectedOutput << 1,0;
    ASSERT_EQ(Hungarian<double>::findMatching(input).indices(),expectedOutput);
}

TEST(AHungarianTest, IntegrationTest_IdenticalPermutation) {
    auto ev = TestMolecules::eightElectrons::square.electrons();

    auto nAlpha = ev.typesVector().countOccurence(Spin::alpha);
    auto nBeta = ev.typesVector().countOccurence(Spin::beta);
    
    // Add random noise
    //double scalingFactor = 0.01;
    //evp.positionsVector().dataRef() += (Eigen::VectorXd::Random(evp.numberOfEntities()*3)*scalingFactor);

    // Create random permutations
    std::random_device rd;
    std::mt19937 g(rd());

    Eigen::PermutationMatrix<Eigen::Dynamic> permAlpha(nAlpha),permBeta(nBeta);
    permAlpha.setIdentity();
    permBeta.setIdentity();

    std::shuffle(permAlpha.indices().data(), permAlpha.indices().data()+permAlpha.indices().size(),g);
    std::shuffle(permBeta.indices().data(), permBeta.indices().data()+permBeta.indices().size(),g);


    auto pAlpha = PositionsVector(ev.positionsVector().asEigenVector().segment(0*3,nAlpha*3));
    auto tAlpha = SpinTypesVector(ev.typesVector().asEigenVector().segment(0,nAlpha));
    auto pBeta = PositionsVector(ev.positionsVector().asEigenVector().segment(nAlpha*3,nBeta*3));
    auto tBeta = SpinTypesVector(ev.typesVector().asEigenVector().segment(nAlpha,nBeta));

    auto evAlpha = ElectronsVector(pAlpha, tAlpha);
    auto evBeta = ElectronsVector(pBeta, tBeta);

    // Permute
    pAlpha.permute(permAlpha);
    tAlpha.permute(permAlpha);
    pBeta.permute(permBeta);
    tBeta.permute(permBeta);

    auto evpAlpha = ElectronsVector(pAlpha, tAlpha);
    auto evpBeta = ElectronsVector(pBeta, tBeta);

    // Find bestmatch permutation to permute ev to evp
    auto costMatrixAlpha = Metrics::positionalDistances(evAlpha.positionsVector(),evpAlpha.positionsVector());
    auto costMatrixBeta = Metrics::positionalDistances(evBeta.positionsVector(),evpBeta.positionsVector());

    auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);
    auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

    // Permute
    pAlpha.permute(bestMatchAlpha.inverse());
    tAlpha.permute(bestMatchAlpha.inverse());
    pBeta.permute(bestMatchBeta.inverse());
    tBeta.permute(bestMatchBeta.inverse());

    Eigen::VectorXd combinedp = Eigen::VectorXd(ev.numberOfEntities()*3);
    Eigen::VectorXi combinedt = Eigen::VectorXi(ev.numberOfEntities());

    combinedp.head(nAlpha*3) = pAlpha.asEigenVector();
    combinedp.tail(nBeta*3) = pBeta.asEigenVector();

    combinedt.head(nAlpha) = tAlpha.asEigenVector();
    combinedt.tail(nBeta) = tBeta.asEigenVector();

    auto evp = ElectronsVector(PositionsVector(combinedp),SpinTypesVector(combinedt));
    ASSERT_EQ(ev,evp);
}
