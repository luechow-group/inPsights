//
// Created by Michael Heuer on 23.05.18.
//
#include <gtest/gtest.h>
#include "StructuralSimilarity.h"
#include "TestMolecules.h"
#include "PositionsVectorTransformer.h"

class AStructuralSimilarityTest : public ::testing::Test {
public:
    double eps = std::numeric_limits<double>::epsilon()*1e3;
    double regularizationParameter = 1.0;

    void SetUp() override {
        ExpansionSettings::defaults();
        ExpansionSettings::mode = ExpansionSettings::Mode::TypeSpecific;
        ParticleKit::create({{Element::H,2},{Element::He,2}},{2,2});
    }
};

TEST_F(AStructuralSimilarityTest , Identity) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_NEAR(StructuralSimilarity::stucturalSimilarity(A,A,regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest , TranslationalSymmetry) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::translated;
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::stucturalSimilarity(A,B,regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry1) {//TODO FIX THIS PROBLEM!!
    auto A = TestMolecules::twoElectrons::oppositeSpin;
    auto B = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ExpansionSettings::Radial::nmax =1;
    ExpansionSettings::Angular::lmax =1;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::stucturalSimilarity(A,B,regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry2) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::permuted2;
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::stucturalSimilarity(A,B,regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest , RotationalSymmetry) {
    auto mol = TestMolecules::H2::ElectronsOutsideCores::offCenter;

    unsigned n = 13;
    
    for (unsigned i = 0; i < n; ++i) {
        double angle = 2*M_PI*double(i)/double(n-1);
        auto pos = mol.electrons().positionsVector();
        PositionsVectorTransformer::rotateAroundAxis(pos,angle,
                                                     mol.atoms()[0].position(),
                                                     mol.atoms()[1].position());

        ElectronsVector rotatedElectrons(pos,mol.electrons().typesVector());
        MolecularGeometry molRotated = {mol.atoms(),rotatedElectrons};

        ASSERT_NEAR(StructuralSimilarity::stucturalSimilarity(mol,molRotated,regularizationParameter),1.0,eps);

    }
}
