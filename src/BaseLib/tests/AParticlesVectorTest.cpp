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
        Particle<Spins::SpinType > e1 = {Spins::SpinType::alpha,{1, 2, 3}};
        Particle<Spins::SpinType > e2 = {Spins::SpinType::alpha,{1, 2, 3}};
        Particle<Spins::SpinType > e3 = {Spins::SpinType::beta ,{4, 5, 6}};
        electrons.append(e1);
        electrons.append(e2);
        electrons.append(e3);

        Particle<Elements::ElementType> a1 = {Elements::ElementType::H ,{1, 2, 3}};
        Particle<Elements::ElementType> a2 = {Elements::ElementType::Og,{4, 5, 6}};
        atoms.append(a1);
        atoms.append(a2);
    };
};

TEST_F(AParticlesVectorTest, BraceInitialization) {

    ParticlesVector<Elements::ElementType> particlesVector(
            {{Elements::ElementType::H , {1, 2, 3}},
             {Elements::ElementType::Og, {4, 5, 6}}}
    );
}

TEST_F(AParticlesVectorTest, CopyConstructor) {
    EXPECT_TRUE(false);
}

TEST_F(AParticlesVectorTest, SpinTypeParticlesVector) {
    std::stringstream stringstream;
    stringstream << electrons;

    std::string expectedOutput = " 1 ea   1.00000   2.00000   3.00000\n"
                                 " 2 ea   1.00000   2.00000   3.00000\n"
                                 " 3 eb   4.00000   5.00000   6.00000\n";
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

TEST_F(AParticlesVectorTest, CountTypeOccurence) {
    ASSERT_EQ(electrons.typesVector().countOccurence(Spins::SpinType::alpha),2);
    ASSERT_EQ(atoms.typesVector().countOccurence(Elements::ElementType::H),1);
}
