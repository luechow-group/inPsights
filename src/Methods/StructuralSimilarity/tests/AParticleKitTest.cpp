//
// Created by Michael Heuer on 08.05.18.
//

#include <gtest/gtest.h>
#include "ParticleKit.h"
#include <iomanip>

using namespace testing;

class AParticleKitTest : public ::testing::Test {
public:
    MolecularGeometry molecularGeometry;
    void SetUp() override {
        Particle<Spins::SpinType > e1 = {Spins::SpinType::alpha,{1, 2, 3}};
        Particle<Spins::SpinType > e2 = {Spins::SpinType::alpha,{1, 2, 3}};
        Particle<Spins::SpinType > e3 = {Spins::SpinType::beta ,{4, 5, 6}};
        molecularGeometry.electrons().append(e1);
        molecularGeometry.electrons().append(e2);
        molecularGeometry.electrons().append(e3);

        Particle<Elements::ElementType> a1 = {Elements::ElementType::He,{1, 2, 3}};
        Particle<Elements::ElementType> a2 = {Elements::ElementType::H ,{4, 5, 6}};
        molecularGeometry.atoms().append(a1);
        molecularGeometry.atoms().append(a2);
    };
};

TEST_F(AParticleKitTest, Constructor1) {
    ParticleKit::create(molecularGeometry.atoms(),0,2);

    // tests order
    ASSERT_EQ(ParticleKit::atomKit[0].first, Elements::ElementType::H);
    ASSERT_EQ(ParticleKit::atomKit[1].first, Elements::ElementType::He);

    ASSERT_EQ(ParticleKit::atomKit[0].second, 1);
    ASSERT_EQ(ParticleKit::atomKit[1].second, 1);

    ASSERT_EQ(ParticleKit::electronKit.first, 2);
    ASSERT_EQ(ParticleKit::electronKit.second, 1);
}

TEST_F(AParticleKitTest, Constructor2) {
    ParticleKit::create(molecularGeometry);

    // tests order
    ASSERT_EQ(ParticleKit::atomKit[0].first, Elements::ElementType::H);
    ASSERT_EQ(ParticleKit::atomKit[1].first, Elements::ElementType::He);

    ASSERT_EQ(ParticleKit::atomKit[0].second, 1);
    ASSERT_EQ(ParticleKit::atomKit[1].second, 1);

    ASSERT_EQ(ParticleKit::electronKit.first, 2);
    ASSERT_EQ(ParticleKit::electronKit.second, 1);
}

TEST_F(AParticleKitTest, WrongMultiplicityOrCharge) {
    EXPECT_DEATH(ParticleKit::create(molecularGeometry.atoms()),"");
}

TEST_F(AParticleKitTest, isSubsetQ) {

    AtomsVector atomsSubset;
    Particle<Elements::ElementType> a1 = {Elements::ElementType::He,{1, 2, 3}};
    atomsSubset.append(a1);

    ParticleKit::create(molecularGeometry.atoms(),0,2);

    ASSERT_TRUE(ParticleKit::isSubsetQ(atomsSubset));
}

TEST_F(AParticleKitTest, NumberedType) {
    ParticleKit::create({{Elements::ElementType::H,2},{Elements::ElementType::Ca,2},{Elements::ElementType::He,2}},{2,3});

    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(0), NumberedElement (Elements::ElementType::H,0));
    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(1), NumberedElement (Elements::ElementType::H,1));
    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(4), NumberedElement (Elements::ElementType::He,0));
    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(5), NumberedElement (Elements::ElementType::He,1));

    ASSERT_EQ(ParticleKit::getNumberedSpinByIndex(1), NumberedSpin(Spins::SpinType::alpha,1));
    ASSERT_EQ(ParticleKit::getNumberedSpinByIndex(4), NumberedSpin(Spins::SpinType::beta,2));
}
