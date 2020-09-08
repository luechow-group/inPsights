// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <Particle.h>
#include <sstream>

using namespace testing;
using namespace Eigen;

class AParticleTest : public Test {
public:
    Eigen::Vector3d pos1{1,2,3};
    Electron electron1{Spin::alpha, pos1};
    Atom atom1{Elements::ElementType::H, pos1};

    void SetUp() override {};
};

TEST_F(AParticleTest, DefaultConstructor) {
    Electron electron2;
    ASSERT_EQ(electron2.position(), Eigen::Vector3d::Zero());
    ASSERT_EQ(electron2.type(), Spins::SpinType::none);
}

TEST_F(AParticleTest, Constructor) {
    Electron electron2(Spin::alpha,pos1);
    ASSERT_EQ(electron1.position(), electron2.position());
}

TEST_F(AParticleTest, CopyConstructor) {
    Electron electron1Copy(electron1);
    ASSERT_EQ(electron1.position(), electron1Copy.position());
}

TEST_F(AParticleTest, Electron) {
    std::stringstream stringstream;
    stringstream << electron1;

    ASSERT_EQ(stringstream.str(), "ea   1.00000   2.00000   3.00000");
    ASSERT_EQ(electron1.type(), Spin::alpha);
}

TEST_F(AParticleTest, Atom) {
    std::stringstream stringstream;
    stringstream << atom1;

    ASSERT_EQ(stringstream.str(), "H    1.00000   2.00000   3.00000");
    ASSERT_EQ(atom1.type(), Elements::ElementType::H);
}

TEST_F(AParticleTest, Charge) {
    ASSERT_EQ(atom1.charge(),+1);
    ASSERT_EQ(electron1.charge(),-1);
}
