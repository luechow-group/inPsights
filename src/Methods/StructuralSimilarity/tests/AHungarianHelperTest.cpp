//
// Created by Michael Heuer on 06.09.18.
//

#include <gmock/gmock.h>
#include <HungarianHelper.h>
#include <TestMolecules.h>
#include <limits>

TEST(AHungarianHelperTest, CombinePermutations){
    Eigen::PermutationMatrix<Eigen::Dynamic>p1,p2;
    p1.setIdentity(4);
    p2.setIdentity(4);

    Eigen::VectorXi expected(8);
    expected << 0,1,2,3,4,5,6,7;

    ASSERT_TRUE(p1.indices().base().isApprox(expected.segment(0,4)));
    ASSERT_TRUE(p2.indices().base().isApprox(expected.segment(0,4)));

    auto p = HungarianHelper::combinePermutations(p1,p2);
    ASSERT_TRUE(p.indices().base().isApprox(expected));
}

TEST(AHungarianHelperTest, SpinSpecificHungarian){
    auto eNormal = TestMolecules::eightElectrons::square.electrons();

    Eigen::VectorXi alphaBetaPositionFlip(8);
    alphaBetaPositionFlip << 4,5,6,7,0,1,2,3;
    Eigen::PermutationMatrix<Eigen::Dynamic> p(alphaBetaPositionFlip);

    auto eFlipped = eNormal;
    eFlipped.positionsVector().permute(p);

    //auto bestMatch = HungarianHelper::spinSpecificHungarian(eNormal, eFlipped, false);
    auto bestMatchFlipped = HungarianHelper::spinSpecificBestMatch(eNormal, eFlipped, true);
    auto bestMatchFlippedInverse = Eigen::PermutationMatrix<Eigen::Dynamic>(bestMatchFlipped.inverse());

    ASSERT_TRUE(bestMatchFlippedInverse.indices().base().isApprox(p.indices().base()));
}

TEST(AHungarianHelperTest, BestMatchNorm) {
    ElectronsVector v1({
        {Spin::alpha, {0,1,2}},
        {Spin::alpha, {0,0,0}}});

    ElectronsVector v2({
            {Spin::alpha, {0,0,0}},
            {Spin::alpha, {0,4,6}}});


    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 1,0;

    auto d2 = Metrics::bestMatchNorm<2,2>(v1, v2);
    ASSERT_EQ(d2.first,5.0);
    ASSERT_TRUE(d2.second.indices().isApprox(expectedPerm));

    auto dInf = Metrics::bestMatchNorm<Eigen::Infinity,2>(v1, v2);
    ASSERT_EQ(dInf.first,5.0);
    ASSERT_TRUE(dInf.second.indices().isApprox(expectedPerm));

    auto d2spin = Metrics::spinSpecifcBestMatchNorm<2,2>(v1, v2);
    ASSERT_EQ(d2spin.first,5.0);
    ASSERT_TRUE(d2spin.second.indices().isApprox(expectedPerm));

    auto dInfInf = Metrics::bestMatchNorm<Eigen::Infinity,Eigen::Infinity>(v1, v2);
    ASSERT_EQ(dInfInf.first,4.0);
    ASSERT_TRUE(dInfInf.second.indices().isApprox(expectedPerm));
}

TEST(AHungarianHelperTest, SpinSpecificBestMatchNorm) {
    ElectronsVector v1({
        {Spin::alpha, {0,1,2}},
        {Spin::beta, {0,0,0}}});

    ElectronsVector v2({
        {Spin::alpha, {0,0,0}},
        {Spin::beta, {0,4,6}}});

    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 0,1;

    auto eps = std::numeric_limits<double>::epsilon()*10;

    auto d2 = Metrics::spinSpecifcBestMatchNorm<2,2>(v1, v2);
    ASSERT_NEAR(d2.first,std::sqrt(1+2*2+4*4+6*6), eps);
    ASSERT_TRUE(d2.second.indices().isApprox(expectedPerm));

    auto dinf = Metrics::spinSpecifcBestMatchNorm<Eigen::Infinity,2>(v1, v2);
    ASSERT_NEAR(dinf.first,std::sqrt(4*4+6*6), eps);
    ASSERT_TRUE(dinf.second.indices().isApprox(expectedPerm));

    auto dInfInf = Metrics::spinSpecifcBestMatchNorm<Eigen::Infinity,Eigen::Infinity>(v1, v2);
    ASSERT_EQ(dInfInf.first,6);
    ASSERT_TRUE(dInfInf.second.indices().isApprox(expectedPerm));
}