//
// Created by Michael Heuer on 08.05.18.
//

#include <gtest/gtest.h>
#include "ParticleKit.h"
#include <iomanip>

using namespace testing;

class AParticlePoolTest : public ::testing::Test {
public:
    MolecularGeometry molecularGeometry;
    void SetUp() override {
        Particle<Spins::SpinType > e1 = {{1,2,3},Spins::SpinType::alpha};
        Particle<Spins::SpinType > e2 = {{1,2,3},Spins::SpinType::alpha};
        Particle<Spins::SpinType > e3 = {{4,5,6},Spins::SpinType::beta};
        molecularGeometry.electrons().append(e1);
        molecularGeometry.electrons().append(e2);
        molecularGeometry.electrons().append(e3);

        Particle<Elements::ElementType> a1 = {{1,2,3},Elements::ElementType::He};
        Particle<Elements::ElementType> a2 = {{4,5,6},Elements::ElementType::H};
        molecularGeometry.atoms().append(a1);
        molecularGeometry.atoms().append(a2);
    };
};

TEST_F(AParticlePoolTest, Constructor1) {
    ParticleKit::create(molecularGeometry.atoms(),0,2);

    // tests order
    ASSERT_EQ(ParticleKit::atomKit[0].first, Elements::ElementType::H);
    ASSERT_EQ(ParticleKit::atomKit[1].first, Elements::ElementType::He);

    ASSERT_EQ(ParticleKit::atomKit[0].second, 1);
    ASSERT_EQ(ParticleKit::atomKit[1].second, 1);

    ASSERT_EQ(ParticleKit::electronKit.first, 2);
    ASSERT_EQ(ParticleKit::electronKit.second, 1);
}

TEST_F(AParticlePoolTest, Constructor2) {
    ParticleKit::create(molecularGeometry);

    // tests order
    ASSERT_EQ(ParticleKit::atomKit[0].first, Elements::ElementType::H);
    ASSERT_EQ(ParticleKit::atomKit[1].first, Elements::ElementType::He);

    ASSERT_EQ(ParticleKit::atomKit[0].second, 1);
    ASSERT_EQ(ParticleKit::atomKit[1].second, 1);

    ASSERT_EQ(ParticleKit::electronKit.first, 2);
    ASSERT_EQ(ParticleKit::electronKit.second, 1);
}

TEST_F(AParticlePoolTest, WrongMultiplicityOrCharge) {
    EXPECT_DEATH(ParticleKit::create(molecularGeometry.atoms()),"");
}

TEST_F(AParticlePoolTest, isSubsetQ) {

    AtomsVector atomsSubset;
    Particle<Elements::ElementType> a1 = {{1,2,3}, Elements::ElementType::He};
    atomsSubset.append(a1);

    ParticleKit::create(molecularGeometry.atoms(),0,2);

    ASSERT_TRUE(ParticleKit::isSubsetQ(atomsSubset));

}