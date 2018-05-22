//
// Created by Michael Heuer on 08.05.18.
//

#include <gtest/gtest.h>
#include "ParticlePool.h"
#include <iomanip>

using namespace testing;

class AParticlePoolTest : public ::testing::Test {
public:

    ElectronsVector electrons;
    AtomsVector atoms;
    void SetUp() override {
        Particle<Spins::SpinType > e1 = {{1,2,3},Spins::SpinType::alpha};
        Particle<Spins::SpinType > e2 = {{1,2,3},Spins::SpinType::alpha};
        Particle<Spins::SpinType > e3 = {{4,5,6},Spins::SpinType::beta};
        electrons.append(e1);
        electrons.append(e2);
        electrons.append(e3);

        Particle<Elements::ElementType> a1 = {{1,2,3},Elements::ElementType::He};
        Particle<Elements::ElementType> a2 = {{4,5,6},Elements::ElementType::H};
        atoms.append(a1);
        atoms.append(a2);
    };
};

TEST_F(AParticlePoolTest, Constructor1) {
    ParticlePool::create(atoms,0,2);

    // tests order
    ASSERT_EQ(ParticlePool::atomKit[0].first, Elements::ElementType::H);
    ASSERT_EQ(ParticlePool::atomKit[1].first, Elements::ElementType::He);

    ASSERT_EQ(ParticlePool::atomKit[0].second, 1);
    ASSERT_EQ(ParticlePool::atomKit[1].second, 1);

    ASSERT_EQ(ParticlePool::electronKit.first, 2);
    ASSERT_EQ(ParticlePool::electronKit.second, 1);
}

TEST_F(AParticlePoolTest, Constructor2) {
    ParticlePool::create(atoms,electrons);

    // tests order
    ASSERT_EQ(ParticlePool::atomKit[0].first, Elements::ElementType::H);
    ASSERT_EQ(ParticlePool::atomKit[1].first, Elements::ElementType::He);

    ASSERT_EQ(ParticlePool::atomKit[0].second, 1);
    ASSERT_EQ(ParticlePool::atomKit[1].second, 1);

    ASSERT_EQ(ParticlePool::electronKit.first, 2);
    ASSERT_EQ(ParticlePool::electronKit.second, 1);
}

TEST_F(AParticlePoolTest, WrongMultiplicityOrCharge) {
    EXPECT_DEATH(ParticlePool::create(atoms),"");
}

TEST_F(AParticlePoolTest, isSubsetQ) {

    AtomsVector atomsSubset;
    Particle<Elements::ElementType> a1 = {{1,2,3},Elements::ElementType::He};
    atomsSubset.append(a1);

    ParticlePool::create(atoms,0,2);

    ASSERT_TRUE(ParticlePool::isSubsetQ(atomsSubset));

}