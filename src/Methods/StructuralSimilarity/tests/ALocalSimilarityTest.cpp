//
// Created by Michael Heuer on 09.05.18.
//


#include <gtest/gtest.h>
#include <ParticlesVector.h>
#include "LocalSimilarity.h"
#include "ParticlePool.h"
#include "ExpansionSettings.h"
#include "Environment.h"

class ALocalSimilarityTest : public ::testing::Test {
public:

    AtomsVector atoms;
    double absError = 1e-8;
    void SetUp() override {
            Atom a0 = {{0,0, 0}, Elements::ElementType::C};
            Atom a1 = {{0,0, 1}, Elements::ElementType::O};
            Atom a2 = {{0,0,-1}, Elements::ElementType::O};
            atoms.append(a0);
            atoms.append(a1);
            atoms.append(a2);
    };
};

TEST_F(ALocalSimilarityTest , NormalizationInGenericMode) {
    ParticlePool pool(atoms);
    ExpansionSettings::defaults();

    Environment e1(atoms,1);
    Environment e2(atoms,2);

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e1,e1),1.0, absError);
    ASSERT_NEAR(LocalSimilarity::localSimilarity(e2,e2),1.0, absError);
};

TEST_F(ALocalSimilarityTest , Diff) {
    ParticlePool pool(atoms);
    ExpansionSettings::defaults();

    Environment e1(atoms,1);
    Environment e2(atoms,2);

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e1,e2),1.0, absError);
};

TEST_F(ALocalSimilarityTest , Cross) {
    ParticlePool pool(atoms);
    ExpansionSettings::defaults();
    ExpansionSettings::Cutoff::cutoffRadius = 1.2;
    ExpansionSettings::Angular::lmax = 3;

    Environment e0(atoms,0);
    Environment e1(atoms,1);
    Environment e2(atoms,2);

    ASSERT_LE(LocalSimilarity::localSimilarity(e0,e1),1.0);
    ASSERT_GE(LocalSimilarity::localSimilarity(e0,e1),0.0);

    ASSERT_LE(LocalSimilarity::localSimilarity(e0,e2),1.0);
    ASSERT_GE(LocalSimilarity::localSimilarity(e0,e2),0.0);

    ASSERT_NEAR(LocalSimilarity::localSimilarity(e0,e1),LocalSimilarity::localSimilarity(e0,e2), absError);
};