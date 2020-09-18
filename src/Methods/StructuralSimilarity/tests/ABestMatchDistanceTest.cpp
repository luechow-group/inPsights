// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <BestMatch.h>
#include <DistanceBasedMetric.h>
#include <TestMolecules.h>
#include <ParticleSelection.h>

TEST(ABestMatchDistanceTest, findTypeSeparatingPermutation) {

    AtomsVector nuclei({{Element::H, {0, 0,-2}},
                        {Element::He,{0, 0,-1}},
                        {Element::Li,{0, 0, 0}},
                        {Element::He,{0, 0, 1}},
                        {Element::H, {0, 0, 2}}});


    auto result = Permutations::findTypeSeparatingPermutation<Element>(nuclei);

    Eigen::VectorXi expected(5);
    expected << 0,4,1,3,2;

    ASSERT_TRUE(result.indices().isApprox(expected));
}

TEST(ABestMatchDistanceTest, TypeSpecificHungarian_Case1) {
    AtomsVector B({
        {Element::H, {0, 0,-3}},
        {Element::He,{0, 0,-2}},
        {Element::Li,{0, 0,-1}},
        {Element::Li,{0, 0, 0}},
        {Element::Li,{0, 0, 1}},
        {Element::He,{0, 0, 2}},
        {Element::H, {0, 0, 3}}});

    AtomsVector A({
        B[6],
        B[5],
        B[4],
        B[3],
        B[2],
        B[1],
        B[0]
    });

    Eigen::VectorXi expected(7);
    expected << 6,5,4,3,2,1,0;

    auto result = Metrics::Similarity::DistanceBased::findTypeSpecificPermutation<Element>(A, B);

    ASSERT_TRUE(result.indices().isApprox(expected));
}

TEST(ABestMatchDistanceTest, TypeSpecificHungarian_Case2) {
    AtomsVector B({
        {Element::Li,{0, 0,-3}},
        {Element::Li,{0, 0,-2}},
        {Element::Li,{0, 0,-1}},
        {Element::He,{0, 0, 0}},
        {Element::He,{0, 0, 1}},
        {Element::He,{0, 0, 2}},
        {Element::H, {0, 0, 3}}});

    AtomsVector A({
        B[3],
        B[2],
        B[0],
        B[1],
        B[4],
        B[5],
        B[6]
    });

    Eigen::VectorXi expected(7);
    expected << 3,2,0,1,4,5,6;

    auto result = Metrics::Similarity::DistanceBased::findTypeSpecificPermutation<Element>(A, B);

    ASSERT_TRUE(result.indices().isApprox(expected));
}

TEST(ABestMatchDistanceTest, TypeSpecificHungarian_Odrered) {
    auto eNormal = TestMolecules::eightElectrons::square.electrons();

    Eigen::VectorXi alphaBetaPositionFlip(8);
    alphaBetaPositionFlip << 4, 5, 6, 7, 0, 1, 2, 3;
    Eigen::PermutationMatrix<Eigen::Dynamic> p(alphaBetaPositionFlip);

    auto eFlipped = eNormal;
    eFlipped.positionsVector().permute(p);
    eFlipped.typesVector().flipSpins();

    auto bestMatchFlipped = Metrics::Similarity::DistanceBased::findTypeSpecificPermutation(eNormal, eFlipped);
    auto bestMatchFlippedInverse = Eigen::PermutationMatrix<Eigen::Dynamic>(bestMatchFlipped.inverse());

    ASSERT_TRUE(bestMatchFlippedInverse.indices().base().isApprox(p.indices().base()));
}

TEST(ABestMatchDistanceTest, BestMatchNorm) {
    ElectronsVector v1({
        {Spin::alpha, {0, 1, 2}},
        {Spin::alpha, {0, 0, 0}}});

    ElectronsVector v2({
        {Spin::alpha, {0, 0, 0}},
        {Spin::alpha, {0, 4, 6}}});


    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 1, 0;

    auto[norm2, perm2] = Metrics::Similarity::DistanceBased::compare<2,2>(
            v1.positionsVector(), v2.positionsVector());
    ASSERT_EQ(norm2, 5.0);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto[normInf, permInf] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
            v1.positionsVector(), v2.positionsVector());
    ASSERT_EQ(normInf, 5.0);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto[normInfInf, permInfInf] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, Eigen::Infinity>(
            v1.positionsVector(), v2.positionsVector());
    ASSERT_EQ(normInfInf, 4.0);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));
}

TEST(ABestMatchDistanceTest, TypeSpecificBestMatchNorm_SameSpin) {
    ElectronsVector v1({
        {Spin::alpha, {0, 1, 2}},
        {Spin::alpha, {0, 0, 0}}});

    ElectronsVector v2({
        {Spin::alpha, {0, 0, 0}},
        {Spin::alpha, {0, 4, 6}}});
    
    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 1, 0;

    auto[norm2, perm2] = Metrics::Similarity::DistanceBased::compare<Spin,2, 2>(v1, v2);
    ASSERT_EQ(norm2, 5.0);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto[normInf, permInf] = Metrics::Similarity::DistanceBased::compare<Spin,Eigen::Infinity, 2>(v1, v2);
    ASSERT_EQ(normInf, 5.0);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto[normInfInf, permInfInf] = Metrics::Similarity::DistanceBased::compare<Spin,Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 4.0);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));
}

TEST(ABestMatchDistanceTest, TypeSpecificBestMatchNorm_DifferentSpin) {
    ElectronsVector v1({
        {Spin::alpha, {0, 1, 2}},
        {Spin::beta,  {0, 0, 0}}});

    ElectronsVector v2({
        {Spin::alpha, {0, 0, 0}},
        {Spin::beta,  {0, 4, 6}}});

    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 0, 1;

    auto eps = std::numeric_limits<double>::epsilon() * 10;

    auto[norm2, perm2] = Metrics::Similarity::DistanceBased::compare<Spin, 2, 2>(v1, v2);
    ASSERT_NEAR(norm2, std::sqrt(1 + 2 * 2 + 4 * 4 + 6 * 6), eps);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto[normInf, permInf] = Metrics::Similarity::DistanceBased::compare<Spin,Eigen::Infinity, 2>(v1, v2);
    ASSERT_NEAR(normInf, std::sqrt(4 * 4 + 6 * 6), eps);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto[normInfInf, permInfInf] = Metrics::Similarity::DistanceBased::compare<Spin,Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 6);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));
}

TEST(ABestMatchDistanceTest, TypeUnspecific_RealMaxima) {
    ElectronsVector v1({
        {Spin::alpha, {-1.924799, -0.000888, -2.199093}},
        {Spin::alpha, {-0.425365, -0.687079, 1.739155}},
        {Spin::alpha, {0.807722,  0.024079,  1.739155}},
        {Spin::alpha, {0.000000,  0.000000,  -1.446226}},
        {Spin::alpha, {0.397920,  -0.688486, -1.772321}},
        {Spin::alpha, {-0.961644, 1.667362,  2.199093}},
        {Spin::alpha, {-0.019931, 0.034495,  -0.660905}},
        {Spin::alpha, {0.000000,  0.000000,  1.446226}},
        {Spin::alpha, {0.961644,  1.667362,  -2.199093}},
        {Spin::beta,  {0.000000,  0.000000,  -1.446226}},
        {Spin::beta,  {0.000000,  0.000000,  1.446226}},
        {Spin::beta,  {-0.963174, -1.666493, 2.199093}},
        {Spin::beta,  {-0.397172, 0.688642,  1.772829}},
        {Spin::beta,  {0.963174,  -1.666493, -2.199093}},
        {Spin::beta,  {-0.807693, -0.025002, -1.739139}},
        {Spin::beta,  {1.924799,  -0.000888, 2.199093}},
        {Spin::beta,  {0.424858,  0.687373,  -1.739140}},
        {Spin::beta,  {0.020247,  -0.035100, 0.660925}}});

    ElectronsVector v2({
        {Spin::alpha, {0.019994,  0.034636,  0.660859}},
        {Spin::alpha, {0.000000,  0.000000,  1.446226}},
        {Spin::alpha, {-0.397902, -0.688437, 1.772419}},
        {Spin::alpha, {0.000000,  0.000000,  -1.446226}},
        {Spin::alpha, {0.961644,  1.667362,  -2.199093}},
        {Spin::alpha, {-0.961644, 1.667362,  2.199093}},
        {Spin::alpha, {0.425527,  -0.687038, -1.739081}},
        {Spin::alpha, {-0.807764, 0.024240,  -1.739086}},
        {Spin::alpha, {1.924799,  -0.000888, 2.199093}},
        {Spin::beta,  {-0.963174, -1.666493, 2.199093}},
        {Spin::beta,  {0.000000,  0.000000,  1.446226}},
        {Spin::beta,  {-0.424878, 0.687387,  1.739079}},
        {Spin::beta,  {0.397266,  0.688792,  -1.772418}},
        {Spin::beta,  {0.000000,  0.000000,  -1.446226}},
        {Spin::beta,  {0.963174,  -1.666493, -2.199093}},
        {Spin::beta,  {0.807718,  -0.025010, 1.739076}},
        {Spin::beta,  {-1.924799, -0.000888, -2.199093}},
        {Spin::beta,  {-0.019967, -0.034674, -0.660818}}});

    auto[norm, perm] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
            v1.positionsVector(),
            v2.positionsVector());

    Eigen::VectorXi expectedPerm(v1.numberOfEntities());
    expectedPerm << 16, 2, 15, 3, 6, 5, 17, 1, 4, 13, 10, 9, 11, 14, 7, 8, 12, 0;

    ASSERT_LT(norm, 0.1);
    ASSERT_TRUE(perm.indices().isApprox(expectedPerm));
}

#include "Metrics.h"

TEST(ABestMatchDistanceTest, NearestElectrons) {
    const MolecularGeometry &BH3 = TestMolecules::BH3::ionic;
    const AtomsVector &nuclei = BH3.atoms();
    const ElectronsVector &electrons = BH3.electrons();
    const ElectronsVector &electrons2 = TestMolecules::BH3::ionicRotated.electrons();
    ElectronsVector electrons3 =
            ElectronsVector({
                                    electrons2[0],
                                    electrons2[1],
                                    electrons2[2],
                                    electrons2[7],
                                    electrons2[4],
                                    electrons2[5],
                                    electrons2[6],
                                    electrons2[3]
                            });

    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.5);

    std::function<double(const Eigen::Vector3d &,
                         const std::vector<Eigen::Vector3d> &)> distanceFunction = Metrics::minimalDistance<2>;
    auto indices1 = ParticleSelection::getNearestElectronsIndices(electrons, nuclei,
                                                                  std::vector<Eigen::Vector3d>({position}), 2, true,
                                                                  100.0,
                                                                  distanceFunction);
    auto indices2 = ParticleSelection::getNearestElectronsIndices(electrons3, nuclei,
                                                                  std::vector<Eigen::Vector3d>({position}), 2, true,
                                                                  100.0,
                                                                  distanceFunction);

    auto[norm, perm] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
            electrons[indices1].positionsVector(),
            electrons3[indices2].positionsVector());

    ASSERT_EQ(norm, 0);
};

TEST(ABestMatchDistanceTest, NearestElectrons2) {
    const MolecularGeometry &BH3 = TestMolecules::BH3::ionic;
    const AtomsVector &nuclei = BH3.atoms();
    const ElectronsVector &electrons = BH3.electrons();
    const ElectronsVector &electrons2 = TestMolecules::BH3::ionicRotated.electrons();

    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 1}, 0.5);

    std::function<double(const Eigen::Vector3d &,
                         const std::vector<Eigen::Vector3d> &)> distanceFunction = Metrics::minimalDistance<2>;
    auto indices1 = ParticleSelection::getNearestElectronsIndices(electrons, nuclei,
                                                                  std::vector<Eigen::Vector3d>({position}), 2, true,
                                                                  100.0,
                                                                  distanceFunction);
    auto indices2 = ParticleSelection::getNearestElectronsIndices(electrons2, nuclei,
                                                                  std::vector<Eigen::Vector3d>({position}), 2, true,
                                                                  100.0,
                                                                  distanceFunction);

    auto[norm, perm] = Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
            electrons[indices1].positionsVector(),
            electrons2[indices2].positionsVector());
    ASSERT_NE(norm, 0);
};
