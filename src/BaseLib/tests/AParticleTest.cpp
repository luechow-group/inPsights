//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include <Particle.h>
#include <sstream>

using namespace testing;
using namespace Eigen;

class AParticleTest : public Test {
public:
    Eigen::Vector3d pos1{1,2,3};
    Electron electron1{pos1,Spins::SpinType::alpha};
    Atom atom1{pos1,Elements::ElementType::H};

    void SetUp() override {};
};

TEST_F(AParticleTest, Constructor) {
    Electron electron2(pos1, Spins::SpinType::alpha);
    ASSERT_EQ(electron1.position(), electron2.position());
}

TEST_F(AParticleTest, CopyConstructor) {
    Electron electron1Copy(electron1);
    ASSERT_EQ(electron1.position(), electron1Copy.position());
}

TEST_F(AParticleTest, SpinTypeParticle) {
    std::stringstream stringstream;
    stringstream << electron1;

    ASSERT_EQ(stringstream.str(), "ea   1.00000   2.00000   3.00000");
    ASSERT_EQ(electron1.type(), Spins::SpinType::alpha);
}

TEST_F(AParticleTest, ElementTypeParticle) {
    std::stringstream stringstream;
    stringstream << atom1;

    ASSERT_EQ(stringstream.str(), "H    1.00000   2.00000   3.00000");
    ASSERT_EQ(atom1.type(), Elements::ElementType::H);
}

/*TEST_F(AParticleTest, Distance) {
    Eigen::Vector3d pos2{1,6,6};
    Electron electron2(pos2, Spins::SpinType::alpha);
    Particle<Elements::ElementType > atom2(pos2, Elements::ElementType::H);

    ASSERT_EQ((Metrics::distance<Spins::SpinType, Spins::SpinType>(electron1,electron2)), 5);
    ASSERT_EQ((Metrics::distance<Elements::ElementType,Elements::ElementType>(atom1,atom2)), 5);
    ASSERT_EQ((Metrics::distance<Elements::ElementType,Spins::SpinType>(atom1,electron2)), 5);
    ASSERT_EQ((Metrics::distance<Elements::ElementType,Spins::SpinType>(atom2,electron2)), 0);
    ASSERT_EQ((Metrics::distance<Elements::ElementType>(atom1,atom2)), 5);
    ASSERT_EQ((Metrics::distance<Spins::SpinType>(electron1,electron2)), 5);
    ASSERT_EQ((Metrics::distance<Elements::ElementType>(atom2,atom2)), 0);

}*/
