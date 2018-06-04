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
        ExpansionSettings::mode = ExpansionSettings::Mode::Chemical;
        ParticleKit::create({{Element::H,2},{Element::He,2}},{2,2});
    }
};

TEST_F(AStructuralSimilarityTest , Identity) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionSettings::Mode::Chemical;
    ParticleKit::create(A);
    
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, A, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest , TranslationalSymmetry) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::translated;
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionSettings::Mode::Chemical;
    ParticleKit::create(A);
    
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry_ReversedOrder) {
    auto A = TestMolecules::twoElectrons::oppositeSpin;
    auto B = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionSettings::Mode::Chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry_FlippedSpins) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::flippedSpins;
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionSettings::Mode::Chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, RotationalSymmetry) {
    auto A = TestMolecules::H2::ElectronsOutsideCores::offCenter;
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionSettings::Mode::Chemical;
    ParticleKit::create(A);

    unsigned n = 13;
    for (unsigned i = 0; i < n; ++i) {
        double angle = 2*M_PI*double(i)/double(n-1);
        auto pos = A.electrons().positionsVector();
        PositionsVectorTransformer::rotateAroundAxis(pos,angle,
                                                     A.atoms()[0].position(),
                                                     A.atoms()[1].position());

        ElectronsVector rotatedElectrons(pos,A.electrons().typesVector());
        MolecularGeometry molRotated = {A.atoms(),rotatedElectrons};

        ASSERT_NEAR(StructuralSimilarity::kernel(A, molRotated, regularizationParameter),1.0,eps);
    }
}

TEST_F(AStructuralSimilarityTest, AlchemicalSimilarity) {
    auto A = TestMolecules::twoElectrons::sameSpinAlpha;
    auto B = TestMolecules::twoElectrons::sameSpinBeta;
    ParticleKit::create({},{2,2});
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));

    ExpansionSettings::mode = ExpansionSettings::Mode::Chemical;
    auto chemical = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_NEAR(chemical, 0.0, eps);

    ExpansionSettings::mode = ExpansionSettings::Mode::Alchemical;
    auto alchemicalSim = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_GT(alchemicalSim, 0.0);
    ASSERT_LT(alchemicalSim, 1.0);
}
