//
// Created by Michael Heuer on 28.08.2018.
//

#include <gmock/gmock.h>
#include <Hungarian.h>
#include <Metrics.h>
#include <TestMolecules.h>

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

#include <algorithm>
#include <random>
TEST(AHungarianTest, IntegrationTest_IdenticalPermutation) {


    //double start = omp_get_wtime();
    auto ev = TestMolecules::eightElectrons::square.electrons();

    auto nAlpha = ev.typesVector().countOccurence(Spin::alpha);
    auto nBeta = ev.typesVector().countOccurence(Spin::beta);
    // assume that vector is ordered
    Interval alphaElectrons({0,nAlpha}), betaElectrons({nAlpha,nBeta});
    // add noise

    auto evp = ev;

    // Add random noise
    double scalingFactor = 0.01;
    //evp.positionsVector().positionsRef() += (Eigen::VectorXd::Random(evp.numberOfEntities()*3)*scalingFactor);

    // Create random permutations
    std::random_device rd;
    std::mt19937 g(rd());

    Eigen::PermutationMatrix<Eigen::Dynamic> permAlpha(nAlpha),permBeta(nBeta);
    permAlpha.setIdentity();
    permBeta.setIdentity();

    std::shuffle(permAlpha.indices().data(), permAlpha.indices().data()+permAlpha.indices().size(),g);
    std::shuffle(permBeta.indices().data(), permBeta.indices().data()+permBeta.indices().size(),g);

    //std::cout << evp << std::endl;
    //std::cout << permAlpha.indices().transpose()<< std::endl;
    //std::cout << permBeta.indices().transpose()<< std::endl;

    // Permute
    evp.slice(alphaElectrons).permute(permAlpha);
    evp.slice(betaElectrons).permute(permBeta);
    //std::cout << evp << std::endl;


    // Find bestmatch permutation to permute ev to evp
    auto costMatrixAlpha = Metrics::positionalDistances(
            PositionsVector(ev.positionsVector().slice(alphaElectrons).dataRef()),
            PositionsVector(evp.positionsVector().slice(alphaElectrons).dataRef()));
    auto costMatrixBeta = Metrics::positionalDistances(
            PositionsVector(ev.positionsVector().slice(betaElectrons).dataRef()),
            PositionsVector(evp.positionsVector().slice(betaElectrons).dataRef()));


    auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);
    auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

    //auto combinedPerm = Eigen::VectorXi(nAlpha+nBeta);
    //combinedPerm.segment(alphaElectrons.start(),alphaElectrons.numberOfEntities()) = bestMatchAlpha.indices();
    //combinedPerm.segment(betaElectrons.start(),betaElectrons.numberOfEntities()) = bestMatchBeta.indices();

    //std::cout << bestMatchAlpha.indices().transpose()<< std::endl;
    //std::cout << bestMatchBeta.indices().transpose()<< std::endl;

    // permute evp to match original ev
    evp.slice(alphaElectrons).permute(bestMatchAlpha.inverse());
    evp.slice(betaElectrons).permute(bestMatchBeta.inverse());

    //std::cout << omp_get_wtime() - start << std::endl;

    //std::cout << evp << std::endl;

    ev.resetSlice(); //TODO WHY DO WE NEED TO DO THIS MANUALLY?
    evp.resetSlice();//TODO WHY DO WE NEED TO DO THIS MANUALLY?

    ASSERT_EQ(ev,evp);

}