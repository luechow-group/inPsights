//
// Created by Michael Heuer on 09.05.18.
//


#include <gtest/gtest.h>
#include <ParticlesVector.h>
#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "ExpansionSettings.h"
#include "Environment.h"
#include "TestMolecules.h"

class ALocalSimilarityTest : public ::testing::Test {
public:

    MolecularGeometry molecule = TestMolecules::CO2::withoutElectrons;
    double absError = std::numeric_limits<double>::epsilon()*1e3;
    void SetUp() override {};
};

TEST_F(ALocalSimilarityTest , GenericNormalization) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e1,e1),1.0, absError);
    ASSERT_NEAR(LocalSimilarity::localSimilarity(e2,e2),1.0, absError);
};

TEST_F(ALocalSimilarityTest , SameEnvironmentsOnDifferentCenters) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e1,e2),1.0, absError);
};

TEST_F(ALocalSimilarityTest , Cross) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();
    ExpansionSettings::Cutoff::cutoffRadius = 1.2;
    ExpansionSettings::Angular::lmax = 3;

    Environment e0(molecule, molecule.atoms()[0].position());
    Environment e1(molecule, molecule.atoms()[1].position());
    Environment e2(molecule, molecule.atoms()[2].position());

    ASSERT_LE(LocalSimilarity::localSimilarity(e0,e1),1.0);
    ASSERT_GE(LocalSimilarity::localSimilarity(e0,e1),0.0);

    ASSERT_LE(LocalSimilarity::localSimilarity(e0,e2),1.0);
    ASSERT_GE(LocalSimilarity::localSimilarity(e0,e2),0.0);

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e0,e1),LocalSimilarity::localSimilarity(e0,e2), absError);
};


TEST_F(ALocalSimilarityTest , TypeSpecificNormalization) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e1,e1),1.0, absError);
    ASSERT_NEAR(LocalSimilarity::localSimilarity(e2,e2),1.0, absError);
};

TEST_F(ALocalSimilarityTest , SameEnvironmentOnDifferentCentersGeneric) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e1,e2),1.0, absError);
};

TEST_F(ALocalSimilarityTest , SameEnvironmentOnDifferentCentersTypeSpecific) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e1,e2),1.0, absError);
};

TEST_F(ALocalSimilarityTest, IsolatedSpecies) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;

    Environment e0(molecule, molecule.atoms()[0].position());
    Environment e1(molecule, molecule.atoms()[1].position());
    Environment e2(molecule, molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e0,e1),0.0, absError);
    ASSERT_NEAR(LocalSimilarity::localSimilarity(e2,e0),0.0, absError);
};

TEST_F(ALocalSimilarityTest , H2) {
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;
    ParticleKit::create(TestMolecules::H2::ElectronsInCores::normal);

    auto H2 = TestMolecules::H2::ElectronsInCores::normal;

    Environment a0(H2, H2.atoms()[0].position());
    Environment e1(H2, H2.electrons()[0].position());

    ASSERT_LE(LocalSimilarity::localSimilarity(a0,e1),1.0);
}
