//
// Created by Michael Heuer on 06.09.18.
//

#include <gmock/gmock.h>
#include <HungarianHelper.h>
#include <TestMolecules.h>

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
    auto bestMatchFlipped = HungarianHelper::spinSpecificHungarian(eNormal, eFlipped, true);
    auto bestMatchFlippedInverse = Eigen::PermutationMatrix<Eigen::Dynamic>(bestMatchFlipped.inverse());

    ASSERT_TRUE(bestMatchFlippedInverse.indices().base().isApprox(p.indices().base()));
    //ASSERT_FALSE(bestMatchFlipped.indices().base().isApprox(p.indices().base()));
    //ASSERT_FALSE(bestMatch.indices().base().isApprox(p.indices().base()));

}

TEST(AHungarianHelperTest, BestMatchNorm) {
    Eigen::VectorXd v1(6),v2(6);
    v1 << 0,1,2,0,0,0;
    v2 << 0,0,0,0,4,6;

    Eigen::VectorXi vp(2);
    vp << 1,0;

    PositionsVector p1(v1);
    PositionsVector p2(v2);
    Eigen::PermutationMatrix<Eigen::Dynamic> perm(vp);

    auto d2 = Metrics::bestMatchNorm<2>(p1,perm,p2);
    ASSERT_EQ(d2,5.0);

    auto dInf = Metrics::bestMatchNorm<Eigen::Infinity>(p1,perm,p2); // Infinty: Pick the largest vector of the vector of positional distances
    ASSERT_EQ(dInf,5.0);
}