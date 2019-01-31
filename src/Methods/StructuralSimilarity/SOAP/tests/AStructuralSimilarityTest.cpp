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
        SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
        ParticleKit::create({{Element::H,2},{Element::He,2}},{2,2});
    }
};

TEST_F(AStructuralSimilarityTest , Identity) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(A);
    
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, A, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest , nmax2) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    Radial::settings.nmax = 2;
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, A, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest , TranslationalSymmetry) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::translated;
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(A);
    
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry_ReversedOrder) {
    auto A = TestMolecules::twoElectrons::oppositeSpin;
    auto B = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry_FlippedSpins) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::flippedSpins;
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, eps);
}

TEST_F(AStructuralSimilarityTest, RotationalSymmetry) {
    auto A = TestMolecules::H2::ElectronsOutsideCores::offCenter;
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(A);

    unsigned n = 13;
    for (unsigned i = 0; i < n; ++i) {
        double angle = 2*M_PI*double(i)/double(n-1);
        auto pos = A.electrons().positionsVector();

        pos.rotateAroundOrigin(angle, A.atoms()[1].position() - A.atoms()[0].position());

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

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    auto chemical = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_NEAR(chemical, 0.0, eps);

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
    auto alchemicalSim = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_GT(alchemicalSim, 0.0);
    ASSERT_LT(alchemicalSim, 1.0);
}

TEST_F(AStructuralSimilarityTest, HeH_H2_Comparison) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::HeH::ElectronsInCores::normal;
    ParticleKit::create({{Element::H,2},{Element::He,1}},{2,2});
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::typeAgnostic;
    auto generic = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_LT(generic, 1.0);
    ASSERT_GT(generic, 0.0);

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    auto chemical = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_LT(chemical, 1.0);
    ASSERT_GT(chemical, 0.0);

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
    auto alchemical = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_LT(alchemical, 1.0);
    ASSERT_GT(alchemical, 0.0);

    // Alchemical similarity must be higher due to the increased similartiy of opposite spins
    ASSERT_LT(chemical,generic);
    ASSERT_LT(chemical,alchemical);
}

TEST_F(AStructuralSimilarityTest, AlchemicalIdentity) {
    auto A = TestMolecules::twoElectrons::sameSpinAlpha;
    auto B = TestMolecules::twoElectrons::sameSpinBeta;
    auto C = TestMolecules::twoElectrons::oppositeSpin;
    ParticleKit::create({},{2,2});
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));

    // Force alchemical identity
    Settings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
    auto ab = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_NEAR(ab, 1.0, eps);

    auto bc = StructuralSimilarity::kernel(B, C, regularizationParameter);
    ASSERT_NEAR(bc, 1.0, eps);

    auto ac = StructuralSimilarity::kernel(A, C, regularizationParameter);
    ASSERT_NEAR(ac, 1.0, eps);
}
