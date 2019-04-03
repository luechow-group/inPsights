//
// Created by heuer on 03.04.19.
//

#include <gmock/gmock.h>
#include <BestMatchSimilarity.h>
#include <ExpansionSettings.h>
#include <TestMolecules.h>

TEST(ABestMatchSimilarityTest, threeElectrons) {
    MolecularGeometry A = {
            AtomsVector(),
            ElectronsVector({{Spin::alpha, {0, 0, -0.5}},
                             {Spin::alpha, {0, 0, 0.5}},
                             {Spin::beta,  {0, 0, 0.3}}
                            })};
    MolecularGeometry B = {
            AtomsVector(),
            ElectronsVector({{Spin::beta,  {0, 0, 0.3}},
                             {Spin::alpha, {0, 0, 0.5}},
                             {Spin::alpha, {0, 0, -0.5}}
                            })};
    MolecularGeometry C = {
            AtomsVector(),
            ElectronsVector({{Spin::alpha, {0, 0, -0.5}},
                             {Spin::alpha, {0, 0, 0.5}},
                             {Spin::beta,  {0, 0, -0.3}}
                            })};
    MolecularGeometry D = {
            AtomsVector(),
            ElectronsVector({{Spin::alpha, {0, 0, -0.5}},
                             {Spin::beta,  {0, 0, -0.3}},
                             {Spin::alpha, {0, 0, 0.5}}
                            })};

    ParticleKit::create(A);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);
    auto specC = MolecularSpectrum(C);
    auto specD = MolecularSpectrum(D);

    auto idPerm = Eigen::PermutationMatrix<Eigen::Dynamic>(A.electrons().numberOfEntities());
    idPerm.setIdentity();


    auto[normAtoA, permAtoA] = BestMatch::Similarity::compare(specA, specA);
    ASSERT_EQ(normAtoA, 1);
    ASSERT_TRUE(permAtoA.indices().isApprox(idPerm.indices()));

    auto[normBtoA, permBtoA] = BestMatch::Similarity::compare(specB, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoA(Eigen::Vector3i(2, 1, 0));
    ASSERT_EQ(normBtoA, 1);
    ASSERT_TRUE(permBtoA.indices().isApprox(refPermBtoA.indices()));

    auto[normCtoA, permCtoA] = BestMatch::Similarity::compare(specC, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoA(Eigen::Vector3i(1, 0, 2));
    ASSERT_EQ(normCtoA, 1);
    ASSERT_TRUE(permCtoA.indices().isApprox(refPermCtoA.indices()));

    auto[normDtoA, permDtoA] = BestMatch::Similarity::compare(specD, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoA(Eigen::Vector3i(1, 2, 0));
    ASSERT_EQ(normDtoA, 1);
    ASSERT_TRUE(permDtoA.indices().isApprox(refPermDtoA.indices()));


    auto[normAtoB, permAtoB] = BestMatch::Similarity::compare(specA, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoB(Eigen::Vector3i(2, 1, 0));
    ASSERT_EQ(normAtoB, 1);
    ASSERT_TRUE(permAtoB.indices().isApprox(refPermAtoB.indices()));

    auto[normBtoB, permBtoB] = BestMatch::Similarity::compare(specB, specB);
    ASSERT_EQ(normBtoB, 1);
    ASSERT_TRUE(permBtoB.indices().isApprox(idPerm.indices()));

    auto[normCtoB, permCtoB] = BestMatch::Similarity::compare(specC, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoB(Eigen::Vector3i(1, 2, 0));
    ASSERT_EQ(normCtoB, 1);
    ASSERT_TRUE(permCtoB.indices().isApprox(refPermCtoB.indices()));

    auto[normDtoB, permDtoB] = BestMatch::Similarity::compare(specD, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoB(Eigen::Vector3i(1, 0, 2));
    ASSERT_EQ(normDtoB, 1);
    ASSERT_TRUE(permDtoB.indices().isApprox(refPermDtoB.indices()));


    auto[normAtoC, permAtoC] = BestMatch::Similarity::compare(specA, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoC(Eigen::Vector3i(1, 0, 2));
    ASSERT_EQ(normAtoC, 1);
    ASSERT_TRUE(permAtoC.indices().isApprox(refPermAtoC.indices()));

    auto[normBtoC, permBtoC] = BestMatch::Similarity::compare(specB, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoC(Eigen::Vector3i(2, 0, 1));
    ASSERT_EQ(normBtoC, 1);
    ASSERT_TRUE(permBtoC.indices().isApprox(refPermBtoC.indices()));

    auto[normCtoC, permCtoC] = BestMatch::Similarity::compare(specC, specC);
    ASSERT_EQ(normCtoC, 1);
    ASSERT_TRUE(permCtoC.indices().isApprox(idPerm.indices()));

    auto[normDtoC, permDtoC] = BestMatch::Similarity::compare(specD, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoC(Eigen::Vector3i(0, 2, 1));
    ASSERT_EQ(normDtoC, 1);
    ASSERT_TRUE(permDtoC.indices().isApprox(refPermDtoC.indices()));


    auto[normAtoD, permAtoD] = BestMatch::Similarity::compare(specA, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoD(Eigen::Vector3i(2, 0, 1));
    ASSERT_EQ(normAtoD, 1);
    ASSERT_TRUE(permAtoD.indices().isApprox(refPermAtoD.indices()));

    auto[normBtoD, permBtoD] = BestMatch::Similarity::compare(specB, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoD(Eigen::Vector3i(1, 0, 2));
    ASSERT_EQ(normBtoD, 1);
    ASSERT_TRUE(permBtoD.indices().isApprox(refPermBtoD.indices()));

    auto[normCtoD, permCtoD] = BestMatch::Similarity::compare(specC, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoD(Eigen::Vector3i(0, 2, 1));
    ASSERT_EQ(normCtoD, 1);
    ASSERT_TRUE(permCtoD.indices().isApprox(refPermCtoD.indices()));

    auto[normDtoD, permDtoD] = BestMatch::Similarity::compare(specD, specD);
    ASSERT_EQ(normDtoD, 1);
    ASSERT_TRUE(permDtoD.indices().isApprox(idPerm.indices()));
}

TEST(ABestMatchSimilarityTest, BH3) {

    auto A = TestMolecules::BH3::ionicMirrored;
    auto B = TestMolecules::BH3::ionic;

    ParticleKit::create(A);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto idPerm = Eigen::PermutationMatrix<Eigen::Dynamic>(A.electrons().numberOfEntities());
    idPerm.setIdentity();

    auto[normAA, permAA] = BestMatch::Similarity::compare(specB, specB);
    ASSERT_EQ(normAA, 1);
    ASSERT_TRUE(permAA.indices().isApprox(idPerm.indices()));

    auto[normBB, permBB] = BestMatch::Similarity::compare(specB, specB);
    ASSERT_EQ(normBB, 1);
    ASSERT_TRUE(permBB.indices().isApprox(idPerm.indices()));

    auto[normAB, permAB] = BestMatch::Similarity::compare(specA, specB);
    Eigen::VectorXi indices(A.electrons().numberOfEntities());
    indices << 0, 1, 4, 5, 2, 3, 6, 7;
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAB(indices);
    ASSERT_EQ(normAB, 1);
    ASSERT_TRUE(permAB.indices().isApprox(refPermAB.indices()));
}