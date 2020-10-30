// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include "StructuralSimilarity.h"
#include "TestMolecules.h"
#include "PositionsVectorTransformer.h"

using namespace SOAP;

class AStructuralSimilarityTest : public ::testing::Test {
public:
    double comparisionEps = SOAP::General::settings.comparisonEpsilon();
    double regularizationParameter = 1.0;

    void SetUp() override {
        General::settings.mode = General::Mode::chemical;
        ParticleKit::create({{Element::H,2},{Element::He,2}},{2,2});
    }
};

TEST_F(AStructuralSimilarityTest , Identity) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(A);
    
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, A, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest , nmax2) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    Radial::settings.nmax = 2;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, A, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest , TranslationalSymmetry) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::translated;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(A);
    
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry_ReversedOrder) {
    auto A = TestMolecules::twoElectrons::oppositeSpin;
    auto B = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest, TShapedGeometry_FlippedSpins) {
    auto A = TestMolecules::fourElectrons::tShaped;
    auto B = TestMolecules::fourElectrons::tShapedSpinFlipped;
    General::settings.mode = General::Mode::chemical;
    General::settings.checkSpinFlip = true;

    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest, PermutationalSymmetry_FlippedSpins) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::H2::ElectronsInCores::flippedSpins;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest, H4linear_Identity) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflected;

    General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest, H6linear_Identity) {
    auto A = TestMolecules::H6::linear::ionicA;
    auto B = TestMolecules::H6::linear::ionicB;

    General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    ParticleKit::create(A);

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0, comparisionEps);
}

TEST_F(AStructuralSimilarityTest, RotationalSymmetry) {
    auto A = TestMolecules::H2::ElectronsOutsideCores::offCenter;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(A);

    unsigned n = 13;
    for (unsigned i = 0; i < n; ++i) {
        double angle = 2*M_PI*double(i)/double(n-1);
        auto pos = A.electrons().positionsVector();

        pos.rotateAroundOrigin(angle, A.atoms()[1].position() - A.atoms()[0].position());

        ElectronsVector rotatedElectrons(pos,A.electrons().typesVector());
        MolecularGeometry molRotated = {A.atoms(),rotatedElectrons};

        ASSERT_NEAR(StructuralSimilarity::kernel(A, molRotated, regularizationParameter),1.0,comparisionEps);
    }
}

TEST_F(AStructuralSimilarityTest, DistortedWater_RotatedAndShiftet) {
    auto A = TestMolecules::Water::distorted::orientation1;
    auto B = TestMolecules::Water::distorted::orientation2;

    General::settings.mode = General::Mode::chemical;
    Radial::settings.nmax = 5;
    Radial::settings.sigmaAtom = 0.25; // sharper basis functions
    Angular::settings.lmax = 5;
    Cutoff::settings.radius = 3.0;
    Cutoff::settings.width = 0.0;

    ParticleKit::create(A);

    double additionalTolerance = 1e-7;  // Reason: geometries themselves only have 4 significant digits
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));
    ASSERT_NEAR(StructuralSimilarity::kernel(A, B, regularizationParameter), 1.0,
                comparisionEps + additionalTolerance);
}

TEST_F(AStructuralSimilarityTest, AlchemicalSimilarity) {
    auto A = TestMolecules::twoElectrons::sameSpinAlpha;
    auto B = TestMolecules::twoElectrons::sameSpinBeta;
    ParticleKit::create({},{2,2});
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));

    General::settings.checkSpinFlip = true;
    General::settings.mode = General::Mode::chemical;
    auto chemicalSpinFlipped = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_EQ(chemicalSpinFlipped, 1.0);

    General::settings.checkSpinFlip = false;
    General::settings.mode = General::Mode::chemical;
    auto chemical = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_EQ(chemical, 0.0);

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 0.5;
    General::settings.mode = General::Mode::alchemical;
    auto alchemicalSim = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_GT(alchemicalSim, 0.0);
    ASSERT_LT(alchemicalSim, 1.0);

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    auto alchemicalSame = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_EQ(alchemicalSame, 1.0);
}

TEST_F(AStructuralSimilarityTest, HeH_H2_Comparison) {
    auto A = TestMolecules::H2::ElectronsInCores::normal;
    auto B = TestMolecules::HeH::ElectronsInCores::normal;
    ParticleKit::create({{Element::H,2},{Element::He,1}},{2,2});
    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));

    General::settings.mode = General::Mode::typeAgnostic;
    auto generic = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_LT(generic, 1.0);
    ASSERT_GT(generic, 0.0);

    General::settings.mode = General::Mode::chemical;
    auto chemical = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_LT(chemical, 1.0);
    ASSERT_GT(chemical, 0.0);

    General::settings.mode = General::Mode::alchemical;
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
    General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;

    General::settings.mode = General::Mode::alchemical;
    auto ab = StructuralSimilarity::kernel(A, B, regularizationParameter);
    ASSERT_NEAR(ab, 1.0, comparisionEps);

    auto bc = StructuralSimilarity::kernel(B, C, regularizationParameter);
    ASSERT_NEAR(bc, 1.0, comparisionEps);

    auto ac = StructuralSimilarity::kernel(A, C, regularizationParameter);
    ASSERT_NEAR(ac, 1.0, comparisionEps);
}
