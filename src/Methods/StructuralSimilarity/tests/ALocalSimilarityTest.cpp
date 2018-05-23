//
// Created by Michael Heuer on 09.05.18.
//


#include <gtest/gtest.h>
#include <ParticlesVector.h>
#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "ExpansionSettings.h"
#include "Environment.h"

class ALocalSimilarityTest : public ::testing::Test {
public:

    MolecularGeometry molecule;
    double absError = 1e-8;
    void SetUp() override {
            Atom a0 = {Elements::ElementType::C,{0,0, 0}};
            Atom a1 = {Elements::ElementType::O,{0,0, 1}};
            Atom a2 = {Elements::ElementType::O,{0,0,-1}};
            molecule.atoms().append(a0);
            molecule.atoms().append(a1);
            molecule.atoms().append(a2);
    };
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

TEST_F(ALocalSimilarityTest , IsolatedSpecies) {
    ParticleKit::create(molecule);
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;

    Environment e0(molecule, molecule.atoms()[0].position());
    Environment e1(molecule, molecule.atoms()[1].position());
    Environment e2(molecule, molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e0,e1),0.0, absError);
    ASSERT_NEAR(LocalSimilarity::localSimilarity(e2,e0),0.0, absError);
};


class ALocalSimilarityTest2 : public ::testing::Test {
public:
    MolecularGeometry molecule;
    void SetUp() override {
        molecule = {
                AtomsVector(
                {{Elements::ElementType::H,{0,0, 0.3705}},
                 {Elements::ElementType::H,{0,0,-0.3705}}}),
                ElectronsVector(
                {{Spins::SpinType::alpha,{0,0, 0.3705}},
                 {Spins::SpinType::beta, {0,0,-0.3705}}})
        };
    }
};

TEST_F(ALocalSimilarityTest2 , Test) {
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;
    ParticleKit::create(molecule);

    Environment a0(molecule, molecule.atoms()[0].position());
    Environment e1(molecule, molecule.electrons()[1].position());

    ASSERT_LE(LocalSimilarity::localSimilarity(a0,e1),1.0);
}
