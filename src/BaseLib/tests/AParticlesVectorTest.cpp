//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include <ParticlesVector.h>
#include <sstream>

using namespace testing;
using namespace Eigen;

class AParticlesVectorTest : public Test {
public:

    ElectronsVector electrons;
    AtomsVector atoms;
    void SetUp() override {
        Particle<Spins::SpinType > e1 = {{1,2,3},Spins::SpinType::alpha};
        Particle<Spins::SpinType > e2 = {{4,5,6},Spins::SpinType::beta};
        electrons.append(e1);
        electrons.append(e2);

        Particle<Elements::ElementType> a1 = {{1,2,3},Elements::ElementType::H};
        Particle<Elements::ElementType> a2 = {{4,5,6},Elements::ElementType::Og};
        atoms.append(a1);
        atoms.append(a2);
    };
};

TEST_F(AParticlesVectorTest, Constructor) {
    EXPECT_TRUE(false);
}

TEST_F(AParticlesVectorTest, CopyConstructor) {
    EXPECT_TRUE(false);
}

TEST_F(AParticlesVectorTest, SpinTypeParticlesVector) {
    std::stringstream stringstream;
    stringstream << electrons;

    std::string expectedOutput = " 1 ea   1.00000   2.00000   3.00000\n"
                                 " 2 eb   4.00000   5.00000   6.00000\n";
    ASSERT_EQ(stringstream.str(), expectedOutput);
}

TEST_F(AParticlesVectorTest, ElementTypeParticlesVector) {
    std::stringstream stringstream;
    stringstream << atoms;

    std::string expectedOutput = " 1 H    1.00000   2.00000   3.00000\n"
                                 " 2 Og   4.00000   5.00000   6.00000\n";
    ASSERT_EQ(stringstream.str(), expectedOutput);
}

TEST_F(AParticlesVectorTest, Distance) {
    EXPECT_TRUE(false);
}
