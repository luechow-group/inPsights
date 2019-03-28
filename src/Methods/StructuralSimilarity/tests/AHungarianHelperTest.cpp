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

    auto [norm2, perm2] = Metrics::bestMatch<2, 2>(v1, v2);
    ASSERT_EQ(norm2, 5.0);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto [normInf, permInf] = Metrics::bestMatch<Eigen::Infinity, 2>(v1, v2);
    ASSERT_EQ(normInf,5.0);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto [normInfInf, permInfInf] = Metrics::bestMatch<Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 4.0);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));

}

TEST(AHungarianHelperTest, SpinSpecificBestMatchNormSameSpin) {
    ElectronsVector v1({
        {Spin::alpha, {0,1,2}},
        {Spin::alpha, {0,0,0}}});

    ElectronsVector v2({
        {Spin::alpha, {0,0,0}},
        {Spin::alpha, {0,4,6}}});


    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 1,0;

    auto [norm2, perm2] = Metrics::spinSpecificBestMatch<2, 2>(v1, v2);
    ASSERT_EQ(norm2, 5.0);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto [normInf, permInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, 2>(v1, v2);
    ASSERT_EQ(normInf,5.0);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto [normInfInf, permInfInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 4.0);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));
}

TEST(AHungarianHelperTest, SpinSpecificBestMatchNormDifferentSpin) {
    ElectronsVector v1({
        {Spin::alpha, {0,1,2}},
        {Spin::beta, {0,0,0}}});

    ElectronsVector v2({
        {Spin::alpha, {0,0,0}},
        {Spin::beta, {0,4,6}}});

    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 0,1;

    auto eps = std::numeric_limits<double>::epsilon()*10;

    auto [norm2, perm2] = Metrics::spinSpecificBestMatch<2, 2>(v1, v2);
    ASSERT_NEAR(norm2, std::sqrt(1+2*2+4*4+6*6), eps);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto [normInf, permInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, 2>(v1, v2);
    ASSERT_NEAR(normInf, std::sqrt(4*4+6*6), eps);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto [normInfInf, permInfInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 6);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));
}

TEST(AHungarianHelperTest, RealMaxima){
    ElectronsVector v1({
    {Spin::alpha,{-1.924799,-0.000888,-2.199093}},
    {Spin::alpha,{-0.425365,-0.687079, 1.739155}},
    {Spin::alpha,{ 0.807722, 0.024079, 1.739155}},
    {Spin::alpha,{ 0.000000, 0.000000,-1.446226}},
    {Spin::alpha,{ 0.397920,-0.688486,-1.772321}},
    {Spin::alpha,{-0.961644, 1.667362, 2.199093}},
    {Spin::alpha,{-0.019931, 0.034495,-0.660905}},
    {Spin::alpha,{ 0.000000, 0.000000, 1.446226}},
    {Spin::alpha,{ 0.961644, 1.667362,-2.199093}},
    {Spin::beta ,{ 0.000000, 0.000000,-1.446226}},
    {Spin::beta ,{ 0.000000, 0.000000, 1.446226}},
    {Spin::beta ,{-0.963174,-1.666493, 2.199093}},
    {Spin::beta ,{-0.397172, 0.688642, 1.772829}},
    {Spin::beta ,{ 0.963174,-1.666493,-2.199093}},
    {Spin::beta ,{-0.807693,-0.025002,-1.739139}},
    {Spin::beta ,{ 1.924799,-0.000888, 2.199093}},
    {Spin::beta ,{ 0.424858, 0.687373,-1.739140}},
    {Spin::beta ,{ 0.020247,-0.035100, 0.660925}}});

    ElectronsVector v2({
    {Spin::alpha,{ 0.019994, 0.034636, 0.660859}},
    {Spin::alpha,{ 0.000000, 0.000000, 1.446226}},
    {Spin::alpha,{-0.397902,-0.688437, 1.772419}},
    {Spin::alpha,{ 0.000000, 0.000000,-1.446226}},
    {Spin::alpha,{ 0.961644, 1.667362,-2.199093}},
    {Spin::alpha,{-0.961644, 1.667362, 2.199093}},
    {Spin::alpha,{ 0.425527,-0.687038,-1.739081}},
    {Spin::alpha,{-0.807764, 0.024240,-1.739086}},
    {Spin::alpha,{ 1.924799,-0.000888, 2.199093}},
    {Spin::beta ,{-0.963174,-1.666493, 2.199093}},
    {Spin::beta ,{ 0.000000, 0.000000, 1.446226}},
    {Spin::beta ,{-0.424878, 0.687387, 1.739079}},
    {Spin::beta ,{ 0.397266, 0.688792,-1.772418}},
    {Spin::beta ,{ 0.000000, 0.000000,-1.446226}},
    {Spin::beta ,{ 0.963174,-1.666493,-2.199093}},
    {Spin::beta ,{ 0.807718,-0.025010, 1.739076}},
    {Spin::beta ,{-1.924799,-0.000888,-2.199093}},
    {Spin::beta ,{-0.019967,-0.034674,-0.660818}}});

    auto [norm, perm] = Metrics::bestMatch<Eigen::Infinity, 2>(v1, v2);

    Eigen::VectorXi expectedPerm(v1.numberOfEntities());
    expectedPerm << 16,2,15,3,6,5,17,1,4,13,10,9,11,14,7,8,12,0;

    ASSERT_LT(norm, 0.1);
    ASSERT_TRUE(perm.indices().isApprox(expectedPerm));
}