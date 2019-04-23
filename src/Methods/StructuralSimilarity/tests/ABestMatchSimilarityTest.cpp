//
// Created by heuer on 03.04.19.
//

#include <gmock/gmock.h>
#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <limits>

using namespace testing;
using namespace SOAP;

class ABestMatchSimilarityTest : public ::testing::Test {
public:
    MolecularGeometry A, B, C, D;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);

        A = {AtomsVector(),
             ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                              {Spin::alpha, {0, 0, 0.5}},
                              {Spin::beta,  {0, 0, 0.3}}})};
        B = {AtomsVector(),
             ElectronsVector({{Spin::beta,  {0, 0, 0.3}},
                              {Spin::alpha, {0, 0, 0.5}},
                              {Spin::alpha, {0, 0, -0.5}}})};
        C = {AtomsVector(),
             ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                              {Spin::alpha, {0, 0, 0.5}},
                              {Spin::beta,  {0, 0,-0.3}}})};
        D = {AtomsVector(),
             ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                              {Spin::beta,  {0, 0,-0.3}},
                              {Spin::alpha, {0, 0, 0.5}}})};
    };
};

TEST_F(ABestMatchSimilarityTest, threeElectrons) {
    ParticleKit::create(A);
    General::settings.mode = General::Mode::chemical;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);
    auto specC = MolecularSpectrum(C);
    auto specD = MolecularSpectrum(D);

    auto idPerm = Eigen::PermutationMatrix<Eigen::Dynamic>(A.electrons().numberOfEntities());
    idPerm.setIdentity();


    auto[normAtoA, permAtoA] = BestMatch::Similarity::compare(specA, specA);
    ASSERT_NEAR(normAtoA, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permAtoA.indices().isApprox(idPerm.indices()));

    auto[normBtoA, permBtoA] = BestMatch::Similarity::compare(specB, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoA(Eigen::Vector3i(2, 1, 0));
    ASSERT_NEAR(normBtoA, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permBtoA.indices().isApprox(refPermBtoA.indices()));

    auto[normCtoA, permCtoA] = BestMatch::Similarity::compare(specC, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoA(Eigen::Vector3i(1, 0, 2));
    ASSERT_NEAR(normCtoA, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permCtoA.indices().isApprox(refPermCtoA.indices()));

    auto[normDtoA, permDtoA] = BestMatch::Similarity::compare(specD, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoA(Eigen::Vector3i(1, 2, 0));
    ASSERT_NEAR(normDtoA, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permDtoA.indices().isApprox(refPermDtoA.indices()));


    auto[normAtoB, permAtoB] = BestMatch::Similarity::compare(specA, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoB(Eigen::Vector3i(2, 1, 0));
    ASSERT_NEAR(normAtoB, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permAtoB.indices().isApprox(refPermAtoB.indices()));

    auto[normBtoB, permBtoB] = BestMatch::Similarity::compare(specB, specB);
    ASSERT_NEAR(normBtoB, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permBtoB.indices().isApprox(idPerm.indices()));

    auto[normCtoB, permCtoB] = BestMatch::Similarity::compare(specC, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoB(Eigen::Vector3i(1, 2, 0));
    ASSERT_NEAR(normCtoB, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permCtoB.indices().isApprox(refPermCtoB.indices()));

    auto[normDtoB, permDtoB] = BestMatch::Similarity::compare(specD, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoB(Eigen::Vector3i(1, 0, 2));
    ASSERT_NEAR(normDtoB, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permDtoB.indices().isApprox(refPermDtoB.indices()));


    auto[normAtoC, permAtoC] = BestMatch::Similarity::compare(specA, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoC(Eigen::Vector3i(1, 0, 2));
    ASSERT_NEAR(normAtoC, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permAtoC.indices().isApprox(refPermAtoC.indices()));

    auto[normBtoC, permBtoC] = BestMatch::Similarity::compare(specB, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoC(Eigen::Vector3i(2, 0, 1));
    ASSERT_NEAR(normBtoC, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permBtoC.indices().isApprox(refPermBtoC.indices()));

    auto[normCtoC, permCtoC] = BestMatch::Similarity::compare(specC, specC);
    ASSERT_NEAR(normCtoC, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permCtoC.indices().isApprox(idPerm.indices()));

    auto[normDtoC, permDtoC] = BestMatch::Similarity::compare(specD, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoC(Eigen::Vector3i(0, 2, 1));
    ASSERT_NEAR(normDtoC, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permDtoC.indices().isApprox(refPermDtoC.indices()));


    auto[normAtoD, permAtoD] = BestMatch::Similarity::compare(specA, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoD(Eigen::Vector3i(2, 0, 1));
    ASSERT_NEAR(normAtoD, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permAtoD.indices().isApprox(refPermAtoD.indices()));

    auto[normBtoD, permBtoD] = BestMatch::Similarity::compare(specB, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoD(Eigen::Vector3i(1, 0, 2));
    ASSERT_NEAR(normBtoD, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permBtoD.indices().isApprox(refPermBtoD.indices()));

    auto[normCtoD, permCtoD] = BestMatch::Similarity::compare(specC, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoD(Eigen::Vector3i(0, 2, 1));
    ASSERT_NEAR(normCtoD, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permCtoD.indices().isApprox(refPermCtoD.indices()));

    auto[normDtoD, permDtoD] = BestMatch::Similarity::compare(specD, specD);
    ASSERT_NEAR(normDtoD, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permDtoD.indices().isApprox(idPerm.indices()));
}

TEST_F(ABestMatchSimilarityTest, BH3) {
    auto a = TestMolecules::BH3::ionicMirrored;
    auto b = TestMolecules::BH3::ionic;

    ParticleKit::create(a);
    General::settings.mode = General::Mode::chemical;

    auto specA = MolecularSpectrum(a);
    auto specB = MolecularSpectrum(b);

    auto idPerm = Eigen::PermutationMatrix<Eigen::Dynamic>(a.electrons().numberOfEntities());
    idPerm.setIdentity();

    auto[normAA, permAA] = BestMatch::Similarity::compare(specB, specB);
    ASSERT_NEAR(normAA, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permAA.indices().isApprox(idPerm.indices()));

    auto[normBB, permBB] = BestMatch::Similarity::compare(specB, specB);
    ASSERT_NEAR(normBB, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permBB.indices().isApprox(idPerm.indices()));

    auto[normAB, permAB] = BestMatch::Similarity::compare(specA, specB);
    Eigen::VectorXi indices(a.electrons().numberOfEntities());
    indices << 0, 1, 4, 5, 2, 3, 6, 7;
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAB(indices);
    ASSERT_NEAR(normAB, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permAB.indices().isApprox(refPermAB.indices()));
}


TEST_F(ABestMatchSimilarityTest, ConvenienceMethods_Unspecific) {
    auto a = TestMolecules::BH3::ionic;
    auto b = TestMolecules::BH3::ionicMirrored;
    b.electrons().typesVector().flipSpins();

    Eigen::VectorXi indices(a.electrons().numberOfEntities());
    indices << 0, 1, 4, 5, 2, 3, 6, 7;
    // careful: many permutations are possible here
    Eigen::PermutationMatrix<Eigen::Dynamic> refPerm(indices);

    auto[norm, perm] = BestMatch::Similarity::compare(b, a, false);
    ASSERT_NEAR(norm, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(perm.indices().isApprox(refPerm.indices()));

    // apply some other permutation that only swaps equivalent electrons of opposite spin
    Eigen::VectorXi otherIndices(a.electrons().numberOfEntities());
    otherIndices << 1, 0, 4, 5, 3, 2, 7, 6;
    Eigen::PermutationMatrix<Eigen::Dynamic> otherPerm(otherIndices);
    b.electrons().typesVector().permute(otherPerm);

    auto[norm2, perm2] = BestMatch::Similarity::compare(b, a, false);
    ASSERT_NEAR(norm2, 1, std::numeric_limits<double>::epsilon());
}

TEST_F(ABestMatchSimilarityTest, ConvenienceMethods_SpinSpecific) {
    auto a = TestMolecules::BH3::ionic;
    auto b = TestMolecules::BH3::ionicMirrored;

    Eigen::VectorXi indices(a.electrons().numberOfEntities());
    indices << 0, 1, 4, 5, 2, 3, 6, 7;

    Eigen::PermutationMatrix<Eigen::Dynamic> refPerm(indices);

    auto[norm, perm] = BestMatch::Similarity::compare(b, a, false);
    ASSERT_NEAR(norm, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(perm.indices().isApprox(refPerm.indices()));
}

TEST_F(ABestMatchSimilarityTest, SpinSpecificWithSpinFlip) {
    auto a = TestMolecules::BH3::ionic;
    auto b = TestMolecules::BH3::ionicMirrored;
    b.electrons().typesVector().flipSpins();

    auto[norm, perm] = BestMatch::Similarity::compare(b, a, true);
    ASSERT_LT(norm, 1);

    Eigen::VectorXi indices(a.electrons().numberOfEntities());
    indices << 0, 1, 4, 5, 2, 3, 6, 7;
    Eigen::PermutationMatrix<Eigen::Dynamic> refPerm(indices);

    auto[normFlipped, permFlipped] = BestMatch::Similarity::compare(b, a, true, true);
    ASSERT_NEAR(normFlipped, 1, std::numeric_limits<double>::epsilon());
    ASSERT_TRUE(permFlipped.indices().isApprox(refPerm.indices()));
}